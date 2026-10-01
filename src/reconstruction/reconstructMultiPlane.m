function [recordedFields, runFolder] = reconstructMultiPlane(inputData, options)
%RECONSTRUCTMULTIPLANE Multi-plane reconstruction with weighted feedback and autofocus.
% inputData contains images, geometry, distanceSteps [m], wavelength [m]
% and pixelSize [m/pixel]. See config/*_example.m for option profiles.
% recordedFields is a cell array of reference-plane complex fields at the
% recorded iterations; sample-plane outputs are saved beneath runFolder.

imgSet = inputData.images;
distanceSteps = inputData.distanceSteps;
settings.wavelength = inputData.wavelength;
settings.pixelSize = inputData.pixelSize;
settings.geometry = inputData.geometry;
settings.iterations = options.iterations;
settings.output = options.output;
settings.runtime.showFigures = options.showFigures;
settings.focus = options.focus;
settings.illumination = sphericalIllumination(size(imgSet{1}), ...
    inputData.geometry, inputData.pixelSize, inputData.wavelength);
recordEvery = options.recordEvery;

nIterations = settings.iterations;
if isempty(recordEvery), recordEvery = nIterations; end
nImages = numel(imgSet);
distanceSteps = distanceSteps(:).';
if numel(distanceSteps) ~= nImages
    error('distanceSteps must contain exactly one value per image.');
end
assert(all(cellfun(@(x) isequal(size(x), size(imgSet{1})), imgSet)), ...
    'All input images must have the same size.');

wavelength = settings.wavelength;
pixelSize = settings.pixelSize;
mnz = settings.geometry;
showFigures = getOption(settings, {'runtime', 'showFigures'}, true);
outputRoot = getOption(settings, {'output', 'rootFolder'}, ...
    fullfile(pwd, 'outputs', 'reconstruction'));
cropOutput = getOption(settings, {'output', 'cropToValidFOV'}, true);
zeroInvalid = getOption(settings, {'output', 'zeroFillInvalid'}, true);
defaultToOriginalFOV = getOption(settings, ...
    {'output', 'defaultToOriginalFOV'}, false);

zPositions = cumsum(distanceSteps);
runStamp = char(datetime('now', 'Format', 'yyyyMMdd_HHmmss_SSS'));
runFolder = fullfile(outputRoot, sprintf('reconstruction_%s_frames%d_iter%d_record%d', ...
    runStamp, nImages, nIterations, recordEvery));
if ~isfolder(runFolder), mkdir(runFolder); end
save(fullfile(runFolder,'run_config.mat'),'options');

hardMasks = getMasks(mnz, 'ValidMaskHard', nImages, size(imgSet{1}), false);
[weightStack, coverageCount] = coverageFusionWeights(hardMasks);
totalWeight = sum(weightStack, 3);
minimumCoverageFraction = getOption(settings, ...
    {'output', 'minimumCoverageFraction'}, 0);
minimumCoverageFraction = min(max(minimumCoverageFraction, 0), 1);
minimumCoverage = getOption(settings, ...
    {'output', 'minimumCoverageCount'}, ...
    max(1, ceil(minimumCoverageFraction * nImages)));
minimumCoverage = min(max(round(minimumCoverage), 1), nImages);
% Rebuild from hard masks, including for existing MAT files with cosine windows.
measurementSupportMask = coverageCount > 0;
outputMask = measurementSupportMask & ...
    coverageCount >= minimumCoverage;
limitHorizontalToReference = getOption(settings, ...
    {'output', 'limitHorizontalToReferenceFOV'}, false);
if limitHorizontalToReference
    referenceColumns = any(hardMasks{1}, 1);
    outputMask = outputMask & repmat(referenceColumns, size(outputMask, 1), 1);
end

initialAmplitude = sqrt(max(imgSet{1}, 0));
initialMask = cast(hardMasks{1}, 'like', initialAmplitude);
background = sum(initialMask(:) .* initialAmplitude(:)) / ...
    max(sum(initialMask(:)), eps('like', initialAmplitude));
initialAmplitude = initialMask .* initialAmplitude + ...
    (1 - initialMask) .* background;
field = initialAmplitude .* settings.illumination;

feedbackA = 0.7;
feedbackB = -0.13 * feedbackA + 0.6;

previous = cell(1, nImages);
previous2 = cell(1, nImages);
rHistory = nan(nIterations, 1);
recordedFields = {};
focusMonitor = [];
focusROI = []; % Selected once per APRW run and reused at recorded iterations.
defaults = defaultReconstructionOptions();
ac = defaults.adaptiveConstraint;
if isfield(options,'adaptiveConstraint')
    names = fieldnames(options.adaptiveConstraint);
    for k = 1:numel(names), ac.(names{k}) = options.adaptiveConstraint.(names{k}); end
end
validateattributes(ac.enabled, {'logical','numeric'}, {'scalar','binary'});
validateattributes(ac.strength, {'numeric'}, {'real','finite','scalar','>=',0,'<=',1});
constraintState = [];
settings.AdaptiveConstraint = ac;

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
        detectorField = propagateAngularSpectrum(field, pixelSize, wavelength, distance);
        predictedAmplitude = abs(detectorField);
        measuredAmplitude = cast(sqrt(max(imgSet{imageIndex}, 0)), ...
            'like', predictedAmplitude);
        measurementMask = cast(hardMasks{imageIndex}, ...
            'like', predictedAmplitude);

        residual = predictedAmplitude - measuredAmplitude;
        errorNumerator = errorNumerator + gatherIfNeeded(sum( ...
            measurementMask(:) .* abs(residual(:)).^2));
        errorDenominator = errorDenominator + gatherIfNeeded(sum( ...
            measurementMask(:) .* abs(measuredAmplitude(:)).^2));

        replacedAmplitude = measurementMask .* measuredAmplitude + ...
            (1 - measurementMask) .* predictedAmplitude;
        detectorField = replacedAmplitude .* exp(1i * angle(detectorField));
        current = propagateAngularSpectrum(detectorField, pixelSize, wavelength, -distance);

        if iteration > 2
            current = (1 + feedbackA + feedbackB) .* current - ...
                feedbackA .* previous{imageIndex} - ...
                feedbackB .* previous2{imageIndex};
        end
        guesses{imageIndex} = current;

        if showFigures && isgraphics(monitor)
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
    weightedField = sum(guessStack .* weightStack, 3) ./ ...
        max(totalWeight, eps('like', totalWeight));

    % Outside the measured support, keep the propagated prediction instead
    % of forcing the padded field to zero. This prevents a hard aperture
    % from producing Fresnel ringing during the final back propagation.
    field = mean(guessStack, 3);
    field(measurementSupportMask) = weightedField(measurementSupportMask);
    constrainedField = min(abs(field), 1.05) .* exp(1i * angle(field));
    field(measurementSupportMask) = constrainedField(measurementSupportMask);

    if ac.enabled && ac.strength > 0
        if isempty(constraintState) && isempty(ac.distance)
            [ac.distance, constraintFocus] = autofocusField( ...
                field,pixelSize,wavelength,settings.focus);
            if strcmp(constraintFocus.roiMode,'manual')
                focusROI = constraintFocus.roiBounds;
            end
            save(fullfile(runFolder,'constraint_autofocus.mat'),'-struct','constraintFocus');
        end
        [field,constraintState,constraintMasks] = applyAdaptiveConstraint( ...
            field,inputData,ac,constraintState);
        ac.whiteAmplitude = constraintState.whiteAmplitude;
        settings.AdaptiveConstraint = ac;
        if iteration == 1
            fprintf('Adaptive constraint: strength %.3g, distance %.6g m, white amplitude %.6g\n', ...
                ac.strength,ac.distance,ac.whiteAmplitude);
        end
    end

    rHistory(iteration) = sqrt(errorNumerator / max(errorDenominator, eps));
    fprintf('Iteration %d/%d: R-factor = %.6f\n', ...
        iteration, nIterations, rHistory(iteration));

    converged = iteration > 2 && ...
        getOption(options, {'stopOnConvergence'}, true) && ...
        all(abs(diff(rHistory(iteration-2:iteration))) < 5e-30);
    shouldRecord = mod(iteration, recordEvery) == 0 || ...
        iteration == nIterations || converged;
    if shouldRecord
        if ac.enabled && ac.strength > 0
            constraintSettings = ac;
            constraintBounds = [constraintState.cols(1),constraintState.rows(1), ...
                constraintState.cols(end),constraintState.rows(end)];
            save(fullfile(runFolder,sprintf('constraint_%04d.mat',iteration)), ...
                'constraintSettings','constraintBounds','constraintMasks');
        end
        recordedFields{end+1, 1} = field; %#ok<AGROW>
        [focusMonitor,focusROI] = saveIteration(field, iteration, runFolder, settings, zPositions, ...
            feedbackA, feedbackB, rHistory, cropOutput, ...
            zeroInvalid, defaultToOriginalFOV, focusMonitor, focusROI);
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
manifest=struct('schemaVersion',1,'algorithm','weighted_multi_plane', ...
    'frameCount',nImages,'iterationsCompleted',iteration, ...
    'latestResult',sprintf('iteration_%04d/results/reconstruction.mat',iteration), ...
    'wavelength_m',wavelength,'pixelSize_m',pixelSize);
fid=fopen(fullfile(runFolder,'run_manifest.json'),'w');
if fid<0,error('Lensless:OutputWrite','Cannot write run manifest.');end
cleanup=onCleanup(@() fclose(fid));
fwrite(fid,jsonencode(manifest,'PrettyPrint',true),'char');clear cleanup;
fprintf('Results saved to: %s\n', runFolder);
end

function [focusMonitor,focusROI] = saveIteration(field, iteration, runFolder, settings, zPositions, ...
        feedbackA, feedbackB, rHistory, ...
        cropOutput, zeroInvalid, ...
        defaultToOriginalFOV, focusMonitor, focusROI)
folder = fullfile(runFolder, sprintf('iteration_%04d', iteration));
resultsFolder = fullfile(folder, 'results');
diagnosticsFolder = fullfile(folder, 'diagnostics');
if ~isfolder(resultsFolder), mkdir(resultsFolder); end
if ~isfolder(diagnosticsFolder), mkdir(diagnosticsFolder); end

if ~isempty(focusROI), settings.focus.roiBounds = focusROI; end
[focusDistance, focusDiagnostics] = autofocusField( ...
    field, settings.pixelSize, settings.wavelength, settings.focus);
if strcmp(focusDiagnostics.roiMode,'manual')
    focusROI = focusDiagnostics.roiBounds;
    settings.focus.roiBounds = focusROI;
end
save(fullfile(diagnosticsFolder, 'autofocus.mat'), '-struct', 'focusDiagnostics');
fprintf('Autofocus distance: %.6g m\n', focusDistance);
objectWithIllumination = propagateAngularSpectrum(field, settings.pixelSize, ...
    settings.wavelength, -focusDistance);

[numRows, numCols] = size(field);
x = ((1:numCols) - floor(numCols/2) - 1) * settings.pixelSize;
y = ((1:numRows) - floor(numRows/2) - 1).' * settings.pixelSize;
[x, y] = meshgrid(x, y);
ledX = settings.geometry.M * settings.pixelSize;
ledY = settings.geometry.N * settings.pixelSize;
ledToSample = settings.geometry.Z - focusDistance;
radius = sqrt((x - ledX).^2 + (y - ledY).^2 + ledToSample^2);
sampleIllumination = exp(1i * (2*pi/settings.wavelength) .* radius);
object = objectWithIllumination .* exp(-1i * angle(sampleIllumination));

geometrySupport = sampleGeometrySupport(size(field), settings.geometry, ...
    zPositions, focusDistance);
sampleCoverageFull = geometrySupport.coverage;
minimumCoverage = getOption(settings, {'output','minimumCoverageCount'}, ...
    max(1,ceil(getOption(settings, {'output','minimumCoverageFraction'},0)*numel(zPositions))));
minimumCoverage = min(max(round(minimumCoverage),1),numel(zPositions));
sampleOutputMask = sampleCoverageFull >= minimumCoverage;
if getOption(settings, {'output','limitHorizontalToReferenceFOV'},false)
    sampleOutputMask = sampleOutputMask & any(geometrySupport.firstMask,1);
end
[safeOutputMask, propagationMargin] = propagationSafeMask( ...
    sampleOutputMask, zPositions, focusDistance, settings);
trustedMask = sampleCoverageFull >= min(2,numel(zPositions)) & safeOutputMask;
[objectOriginalFOV, originalValidMask, originalBounds] = applyOutputMask( ...
    object, geometrySupport.firstMask, true, zeroInvalid);
[object, validMask, bounds] = applyOutputMask( ...
    object, safeOutputMask, cropOutput, zeroInvalid);
trustedMask = trustedMask(bounds(2):bounds(4), bounds(1):bounds(3));
sampleCoverage = sampleCoverageFull(bounds(2):bounds(4), bounds(1):bounds(3));
unionAmplitude = normalizePercentile(abs(object), 1, 99, validMask);
originalAmplitude = normalizePercentile(abs(objectOriginalFOV), 1, 99);
[unionPhaseRGB, phaseDisplayLimits] = phaseHeatmap(angle(object), validMask);
[originalPhaseRGB, originalPhaseDisplayLimits] = ...
    phaseHeatmap(angle(objectOriginalFOV));

if defaultToOriginalFOV
    amplitude = originalAmplitude;
    phaseRGB = originalPhaseRGB;
    defaultView = 'first_camera_sample_fov';
    displayMask = originalValidMask;
    focusDiagnostics.displayBounds = originalBounds;
else
    amplitude = unionAmplitude;
    phaseRGB = unionPhaseRGB;
    defaultView = 'geometric_sample_union';
    displayMask = validMask;
    focusDiagnostics.displayBounds = bounds;
end
save(fullfile(diagnosticsFolder, 'autofocus.mat'), '-struct', 'focusDiagnostics');

save(fullfile(resultsFolder, 'reconstruction.mat'), ...
    'field', 'object', 'objectOriginalFOV', 'validMask', 'trustedMask', ...
    'bounds', 'originalBounds', 'focusDistance', ...
    'propagationMargin', 'phaseDisplayLimits', 'originalPhaseDisplayLimits', ...
    'defaultView', 'sampleCoverage', 'originalValidMask');
save(fullfile(diagnosticsFolder, 'sample_geometry.mat'), 'geometrySupport');
save(fullfile(diagnosticsFolder, 'meta.mat'), 'settings', 'zPositions', ...
    'feedbackA', 'feedbackB', 'rHistory', 'focusDistance', 'bounds', ...
    'originalBounds', 'propagationMargin');

writeMaskedPNG(amplitude, fullfile(resultsFolder, 'amplitude.png'), displayMask);
writeMaskedPNG(phaseRGB, fullfile(resultsFolder, 'phase_heatmap.png'), displayMask);
if ~cropOutput
    % Color fusion needs the full numerical canvas, but each channel's saved
    % iteration also provides standalone grayscale views without padding.
    writeCroppedReconstruction(object,validMask,resultsFolder);
end
writeMaskedPNG(unionAmplitude, fullfile(resultsFolder, 'amplitude_union.png'), validMask);
writeMaskedPNG(unionPhaseRGB, fullfile(resultsFolder, 'phase_heatmap_union.png'), validMask);
imwrite(originalAmplitude, fullfile(resultsFolder, 'originalFOV_amplitude.png'));
imwrite(originalPhaseRGB, ...
    fullfile(resultsFolder, 'originalFOV_phase_heatmap.png'));
% Standard runs update one persistent window. Headless runs export through
% a temporary invisible figure, leaving no windows behind.
showFocus = getOption(settings, {'runtime','showFigures'}, true);
saveFocus = getOption(settings, {'output','saveFocusPlot'}, true);
if showFocus || saveFocus
    visible = 'off';
    if showFocus, visible = 'on'; end
    focusMonitor = plotAutofocus(focusDiagnostics, amplitude, phaseRGB, ...
        sprintf('Iteration %d | reconstructed sample field', iteration), ...
        focusMonitor, visible);
    if ~showFocus, cleanupFigure = onCleanup(@() close(focusMonitor)); end
    if saveFocus
        exportgraphics(focusMonitor, fullfile(diagnosticsFolder, 'autofocus.png'), ...
            'Resolution', 150);
    end
    if ~showFocus
        clear cleanupFigure;
        focusMonitor = [];
    end
end
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

function [safeMask, margin] = propagationSafeMask( ...
        mask, zPositions, focusDistance, settings)
trimBoundary = getOption(settings, ...
    {'output', 'trimPropagationBoundary'}, false);
safeMask = mask;
margin = 0;
if ~trimBoundary || all(mask(:))
    return;
end

% Same finite-sensor support estimate used by padded angular-spectrum models.
maximumDistance = max([abs(zPositions(:)); abs(focusDistance)]);
requestedMargin = ceil(maximumDistance * settings.wavelength / ...
    (2 * settings.pixelSize^2));
distanceToInvalid = bwdist(~mask);
maximumAvailable = max(distanceToInvalid(mask), [], 'all');
margin = min(requestedMargin, max(floor(maximumAvailable) - 1, 0));
if margin > 0
    safeMask = mask & distanceToInvalid > margin;
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
