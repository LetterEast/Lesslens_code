function recordedFields = APRW(inputData, options)
%APRW Multi-plane reconstruction with weighted feedback and autofocus.

imgSet = inputData.images;
distanceSteps = inputData.distanceSteps;
mainPara.WaveLength = inputData.wavelength;
mainPara.PixelSize = inputData.pixelSize;
mainPara.MNZ_result = inputData.geometry;
mainPara.nIterative = options.iterations;
mainPara.ParaD_Sample2CCD.D_Sample2CCDPre = options.focus.prior;
mainPara.ParaD_Sample2CCD.D_Sample2CCDHalfRange = options.focus.halfRange;
mainPara.ParaD_Sample2CCD.rough = options.focus.step;
mainPara.Output = options.output;
mainPara.Runtime.showFigures = options.showFigures;
mainPara.IllumSet = sphericalIllumination(size(imgSet{1}), ...
    inputData.geometry, inputData.pixelSize, inputData.wavelength);
recordEvery = options.recordEvery;

nIterations = mainPara.nIterative;
if nargin < 4 || isempty(recordEvery), recordEvery = nIterations; end
nImages = numel(imgSet);
distanceSteps = distanceSteps(:).';
if numel(distanceSteps) ~= nImages
    error('distanceSteps must contain exactly one value per image.');
end
assert(all(cellfun(@(x) isequal(size(x), size(imgSet{1})), imgSet)), ...
    'All input images must have the same size.');

wavelength = mainPara.WaveLength;
pixelSize = mainPara.PixelSize;
mnz = mainPara.MNZ_result;
showFigures = getOption(mainPara, {'Runtime', 'showFigures'}, true);
outputRoot = getOption(mainPara, {'Output', 'rootFolder'}, ...
    fullfile(pwd, 'ResultFolder'));
cropOutput = getOption(mainPara, {'Output', 'cropToValidFOV'}, true);
zeroInvalid = getOption(mainPara, {'Output', 'zeroFillInvalid'}, true);

zPositions = cumsum(distanceSteps);
runStamp = char(datetime('now', 'Format', 'yyyyMMdd_HHmmss'));
runFolder = fullfile(outputRoot, sprintf('APRW_%s_N%d_I%d_rec%d', ...
    runStamp, nImages, nIterations, recordEvery));
if ~isfolder(runFolder), mkdir(runFolder); end

hardMasks = getMasks(mnz, 'ValidMaskHard', nImages, size(imgSet{1}), false);
softMasks = getMasks(mnz, 'ValidMask', nImages, size(imgSet{1}), true);
outputMask = any(cat(3, hardMasks{:}), 3);
trustedMask = sum(cat(3, hardMasks{:}), 3) >= min(2, nImages);

initialAmplitude = sqrt(max(imgSet{1}, 0));
initialMask = cast(hardMasks{1}, 'like', initialAmplitude);
background = sum(initialMask(:) .* initialAmplitude(:)) / ...
    max(sum(initialMask(:)), eps('like', initialAmplitude));
initialAmplitude = initialMask .* initialAmplitude + ...
    (1 - initialMask) .* background;
field = initialAmplitude .* mainPara.IllumSet;

feedbackA = 0.7;
feedbackB = -0.13 * feedbackA + 0.6;
previous = cell(1, nImages);
previous2 = cell(1, nImages);
rHistory = nan(nIterations, 1);
recordedFields = {};

if showFigures
    monitor = figure('Name', 'APRW reconstruction');
    amplitudeAxes = subplot(1, 2, 1, 'Parent', monitor);
    phaseAxes = subplot(1, 2, 2, 'Parent', monitor);
    amplitudeImage = imshow(zeros(size(field)), [], 'Parent', amplitudeAxes);
    phaseImage = imshow(zeros(size(field)), [], 'Parent', phaseAxes);
end

tic;
for iteration = 1:nIterations
    guesses = cell(1, nImages);
    errorNumerator = 0;
    errorDenominator = 0;

    for imageIndex = 1:nImages
        distance = zPositions(imageIndex);
        detectorField = propagate(field, pixelSize, wavelength, distance);
        predictedAmplitude = abs(detectorField);
        measuredAmplitude = cast(sqrt(max(imgSet{imageIndex}, 0)), ...
            'like', predictedAmplitude);
        measurementMask = cast(hardMasks{imageIndex}, ...
            'like', predictedAmplitude);

        residual = predictedAmplitude - measuredAmplitude;
        errorNumerator = errorNumerator + gather(sum( ...
            measurementMask(:) .* abs(residual(:)).^2));
        errorDenominator = errorDenominator + gather(sum( ...
            measurementMask(:) .* abs(measuredAmplitude(:)).^2));

        replacedAmplitude = measurementMask .* measuredAmplitude + ...
            (1 - measurementMask) .* predictedAmplitude;
        detectorField = replacedAmplitude .* exp(1i * angle(detectorField));
        current = propagate(detectorField, pixelSize, wavelength, -distance);

        if iteration > 2
            current = (1 + feedbackA + feedbackB) .* current - ...
                feedbackA .* previous{imageIndex} - ...
                feedbackB .* previous2{imageIndex};
        end
        guesses{imageIndex} = current;

        if showFigures
            set(amplitudeImage, 'CData', mat2gray(abs(current)));
            set(phaseImage, 'CData', mat2gray(angle(current)));
            title(amplitudeAxes, sprintf('Iteration %d/%d, image %d/%d', ...
                iteration, nIterations, imageIndex, nImages));
            title(phaseAxes, 'Phase');
            drawnow limitrate;
        end
    end

    previous2 = previous;
    previous = guesses;
    guessStack = cat(3, guesses{:});
    weightStack = cat(3, softMasks{:});
    totalWeight = sum(weightStack, 3);
    field = sum(guessStack .* weightStack, 3) ./ ...
        max(totalWeight, eps('like', totalWeight));
    field = min(abs(field), 1.05) .* exp(1i * angle(field));

    rHistory(iteration) = sqrt(errorNumerator / max(errorDenominator, eps));
    fprintf('Iteration %d/%d: R-factor = %.6f\n', ...
        iteration, nIterations, rHistory(iteration));

    converged = iteration > 2 && ...
        all(abs(diff(rHistory(iteration-2:iteration))) < 5e-3);
    shouldRecord = mod(iteration, recordEvery) == 0 || ...
        iteration == nIterations || converged;
    if shouldRecord
        recordedFields{end+1, 1} = field; %#ok<AGROW>
        saveIteration(field, iteration, runFolder, mainPara, zPositions, ...
            feedbackA, feedbackB, rHistory, outputMask, trustedMask, ...
            cropOutput, zeroInvalid);
    end
    if converged
        fprintf('R-factor has converged; stopping after iteration %d.\n', iteration);
        break;
    end
end

if showFigures && isgraphics(monitor)
    figure('Name', 'R-factor convergence');
    plot(1:iteration, rHistory(1:iteration), '-o', 'LineWidth', 1.5);
    xlabel('Iteration'); ylabel('R-factor'); grid on;
end
toc;
fprintf('Results saved to: %s\n', runFolder);
end

function saveIteration(field, iteration, runFolder, mainPara, zPositions, ...
        feedbackA, feedbackB, rHistory, outputMask, trustedMask, ...
        cropOutput, zeroInvalid)
folder = fullfile(runFolder, sprintf('Iter_%04d', iteration));
if ~isfolder(folder), mkdir(folder); end

focusDistance = autofocus(field, mainPara, folder);
objectWithIllumination = propagate(field, mainPara.PixelSize, ...
    mainPara.WaveLength, -focusDistance);

[numRows, numCols] = size(field);
x = ((1:numCols) - floor(numCols/2) - 1) * mainPara.PixelSize;
y = ((1:numRows) - floor(numRows/2) - 1).' * mainPara.PixelSize;
[x, y] = meshgrid(x, y);
ledX = mainPara.MNZ_result.M * mainPara.PixelSize;
ledY = mainPara.MNZ_result.N * mainPara.PixelSize;
ledToSample = mainPara.MNZ_result.Z - focusDistance;
radius = sqrt((x - ledX).^2 + (y - ledY).^2 + ledToSample^2);
sampleIllumination = exp(1i * (2*pi/mainPara.WaveLength) .* radius);
object = objectWithIllumination .* exp(-1i * angle(sampleIllumination));

% Save a fixed-size result at the exact location of the first input image.
% This excludes the low-coverage synthetic field added by later shifts.
[objectOriginalFOV, originalBounds] = cropFirstImageFOV( ...
    object, mainPara.MNZ_result);
save(fullfile(folder, sprintf('Object_originalFOV_iter%04d.mat', iteration)), ...
    'objectOriginalFOV', 'originalBounds');
originalAmplitude = normalizePercentile(abs(objectOriginalFOV), 1, 99);
originalPhase = normalizePercentile(angle(objectOriginalFOV), 1, 99);
imwrite(originalAmplitude, fullfile(folder, sprintf( ...
    'Object_originalFOV_amp_iter%04d.png', iteration)));
imwrite(originalPhase, fullfile(folder, sprintf( ...
    'Object_originalFOV_phs_iter%04d.png', iteration)));

[object, validMask, bounds] = applyOutputMask( ...
    object, outputMask, cropOutput, zeroInvalid);
trustedMask = trustedMask(bounds(2):bounds(4), bounds(1):bounds(3));
trustedObject = object;
trustedObject(~trustedMask) = 0;

save(fullfile(folder, 'ImgRec.mat'), 'field');
save(fullfile(folder, sprintf('Object_iter%04d.mat', iteration)), 'object');
save(fullfile(folder, sprintf('Object_trusted_iter%04d.mat', iteration)), ...
    'trustedObject', 'trustedMask');
save(fullfile(folder, 'meta.mat'), 'mainPara', 'zPositions', ...
    'feedbackA', 'feedbackB', 'rHistory', 'focusDistance', 'bounds', ...
    'originalBounds');
save(fullfile(folder, 'OutputFOV.mat'), 'validMask', 'trustedMask', 'bounds');

amplitude = normalizePercentile(abs(object), 1, 99);
phase = normalizePercentile(angle(object), 1, 99);
imwrite(amplitude, fullfile(folder, sprintf('Object_amp_iter%04d.png', iteration)));
imwrite(phase, fullfile(folder, sprintf('Object_phs_iter%04d.png', iteration)));
end

function [field, bounds] = cropFirstImageFOV(field, mnz)
if ~isfield(mnz, 'padSize') || ~isfield(mnz, 'orig_size')
    bounds = [1, 1, size(field, 2), size(field, 1)];
    return;
end

padY = mnz.padSize(1);
padX = mnz.padSize(2);
originalHeight = mnz.orig_size(1);
originalWidth = mnz.orig_size(2);
firstRow = padY + 1;
firstColumn = padX + 1;
lastRow = firstRow + originalHeight - 1;
lastColumn = firstColumn + originalWidth - 1;
if lastRow > size(field, 1) || lastColumn > size(field, 2)
    error('Stored original-FOV bounds exceed the reconstruction canvas.');
end
bounds = [firstColumn, firstRow, lastColumn, lastRow];
field = field(firstRow:lastRow, firstColumn:lastColumn);
end

function illumination = sphericalIllumination(imageSize, geometry, pixelSize, wavelength)
height = imageSize(1); width = imageSize(2);
x = ((1:width) - floor(width/2) - 1) * pixelSize;
y = ((1:height) - floor(height/2) - 1).' * pixelSize;
[x, y] = meshgrid(x, y);
radius = sqrt((x - geometry.M*pixelSize).^2 + ...
    (y - geometry.N*pixelSize).^2 + geometry.Z^2);
illumination = exp(1i * (2*pi/wavelength) .* radius);
end

function bestDistance = autofocus(field, mainPara, folder)
settings = mainPara.ParaD_Sample2CCD;
distances = settings.D_Sample2CCDPre + ...
    (-settings.D_Sample2CCDHalfRange:settings.rough: ...
    settings.D_Sample2CCDHalfRange);
if isempty(distances), distances = settings.D_Sample2CCDPre; end
focusMetric = zeros(size(distances));
for index = 1:numel(distances)
    object = propagate(field, mainPara.PixelSize, ...
        mainPara.WaveLength, -distances(index));
    amplitude = abs(object);
    highPass = amplitude - imgaussfilt(amplitude, 2);
    focusMetric(index) = sum(abs(fft2(highPass)), 'all');
end
[~, bestIndex] = max(focusMetric);
bestDistance = distances(bestIndex);
save(fullfile(folder, 'autofocus.mat'), ...
    'distances', 'focusMetric', 'bestDistance');
fprintf('Autofocus distance: %.6g m\n', bestDistance);
end

function masks = getMasks(mnz, fieldName, count, imageSize, useSoft)
if isfield(mnz, fieldName) && numel(mnz.(fieldName)) >= count
    masks = mnz.(fieldName)(1:count);
else
    masks = repmat({ones(imageSize)}, 1, count);
end
for index = 1:count
    if useSoft
        masks{index} = min(max(double(masks{index}), 0), 1);
    else
        masks{index} = logical(masks{index});
    end
end
end

function [field, mask, bounds] = applyOutputMask(field, mask, crop, zeroInvalid)
if ~any(mask(:)), error('The output valid-field mask is empty.'); end
if crop
    [rows, cols] = find(mask);
    bounds = [min(cols), min(rows), max(cols), max(rows)];
    field = field(bounds(2):bounds(4), bounds(1):bounds(3));
    mask = mask(bounds(2):bounds(4), bounds(1):bounds(3));
else
    bounds = [1, 1, size(field, 2), size(field, 1)];
end
if zeroInvalid, field(~mask) = 0; end
end

function image = normalizePercentile(data, lowPercentile, highPercentile)
data = double(gather(data));
limits = prctile(data(:), [lowPercentile, highPercentile]);
image = (data - limits(1)) / max(limits(2) - limits(1), eps);
image = min(max(image, 0), 1);
end

function value = getOption(structure, path, defaultValue)
value = structure;
for index = 1:numel(path)
    if ~isstruct(value) || ~isfield(value, path{index})
        value = defaultValue;
        return;
    end
    value = value.(path{index});
end
end
