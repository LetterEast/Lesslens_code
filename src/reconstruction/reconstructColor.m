function outputRoot = reconstructColor(config)
%RECONSTRUCTCOLOR Reconstruct prepared channel MATs and fuse registered RGB.
% config.inputFile contains colorInput.channelFiles (three paths, B/G/R by
% default). Each channel MAT contains the documented APRW inputData struct.
% Channels are loaded one at a time to keep memory bounded.
if ~isfile(config.inputFile)
    error('Lensless:MissingColorInput', ...
        'Prepared color manifest is missing: %s. See docs/INPUT_FORMAT.md.', config.inputFile);
end
manifest = load(config.inputFile, 'colorInput');
if ~isfield(manifest, 'colorInput') || ~isfield(manifest.colorInput, 'channelFiles')
    error('Lensless:InvalidColorInput', 'Manifest must contain colorInput.channelFiles.');
end
channelFiles = manifest.colorInput.channelFiles;
assert(iscell(channelFiles) && numel(channelFiles)==3, 'Exactly three channel MAT paths are required.');
manifestFolder = fileparts(config.inputFile);
if isempty(manifestFolder), manifestFolder = pwd; end
outputRoot = config.outputRoot;
reuseExistingReconstructions = config.reuseExistingReconstructions;
% A new interactive selection must not be bypassed by a cached channel run.
if isfield(config.options.focus,'roiMode') && strcmp(config.options.focus.roiMode,'manual') && ...
        (~isfield(config.options.focus,'roiBounds') || isempty(config.options.focus.roiBounds))
    reuseExistingReconstructions = false;
end
channelNames = config.channelNames;
channelFolders = channelNames; % Retain historical output metadata field.
rgbOrder = config.rgbOrder;
assert(numel(channelNames)==3 && isequal(sort(rgbOrder(:).'), [1,2,3]), ...
    'Three names and an RGB permutation of 1:3 are required.');
assert(iscellstr(channelNames) && numel(unique(channelNames))==3 && ...
    all(~cellfun(@isempty, regexp(channelNames, '^[A-Za-z0-9_-]+$', 'once'))), ...
    'Channel names must be unique simple directory names.');
colorBalance = config.colorBalance;
physicalIntensity = config.physicalIntensity;
inputBrightness = config.inputBrightness;
options = config.options;
grayscaleOutput = options.output;
% Full canvas is needed for physical-coordinate RGB registration. Export
% separate grayscale views using the user's original output settings.
options.output.cropToValidFOV = false;
options.output.zeroFillInvalid = false;
options.inputBrightness = inputBrightness;
if ~isfolder(outputRoot), mkdir(outputRoot); end
channelCount = numel(channelFiles);
channelAmplitude = cell(1, channelCount);
channelMask = cell(1, channelCount);
channelCoverage = cell(1, channelCount);
channelGeometry = cell(1, channelCount);
focusDistances = zeros(1, channelCount);
wavelengths = zeros(1, channelCount);
channelRunFolders = cell(1, channelCount);
channelImageCounts = zeros(1, channelCount);
distanceStepsByChannel = cell(1, channelCount);
inputBrightnessGains = cell(1, channelCount);
inputBrightnessLevels = cell(1, channelCount);

for channel = 1:channelCount
    inputFile = char(channelFiles{channel});
    if isempty(regexp(inputFile, '^([A-Za-z]:[\\/]|[\\/])', 'once'))
        inputFile = fullfile(manifestFolder, inputFile);
    end
    if ~isfile(inputFile)
        error('Lensless:MissingInput', 'Prepared channel MAT is missing: %s', inputFile);
    end
    fprintf('\n[%d/%d] Reconstructing %s channel\n', channel, channelCount, channelNames{channel});
    channelOutputRoot = fullfile(outputRoot, ['channel_' channelNames{channel}]);
    options.output.rootFolder = channelOutputRoot;
    sourceInfo = reconstructionSourceInfo(inputFile, options);
    runFolder = '';
    if reuseExistingReconstructions
        runFolder = latestCompatibleRunFolder(channelOutputRoot, sourceInfo);
    end
    if isempty(runFolder)
        prepared = load(inputFile, 'inputData');
        if ~isfield(prepared, 'inputData')
            error('Lensless:InvalidInput', 'Channel MAT must contain inputData: %s', inputFile);
        end
        inputData = prepared.inputData;
        validateReconstructionInput(inputData);
        frameBrightnessGains = [];
        frameBrightnessLevels = [];
        if inputBrightness.enabled
            [inputData.images, frameBrightnessGains, frameBrightnessLevels] = ...
                equalizeMeasurementBrightness(inputData.images, ...
                inputData.geometry.ValidMaskHard, inputBrightness);
        end
        [~, runFolder] = reconstructMultiPlane(inputData, options);
        save(fullfile(runFolder, 'color_source.mat'), 'sourceInfo', ...
            'frameBrightnessGains', 'frameBrightnessLevels');
        clear prepared inputData;
    else
        fprintf('Reusing matching reconstruction: %s\n', runFolder);
    end
    sourceData = load(fullfile(runFolder, 'color_source.mat'));
    inputBrightnessGains{channel} = sourceData.frameBrightnessGains;
    inputBrightnessLevels{channel} = sourceData.frameBrightnessLevels;
    [resultFile, iterationFolder] = latestReconstructionFile(runFolder);
    metadata = load(fullfile(iterationFolder, 'diagnostics', 'meta.mat'));
    if isfield(metadata,'settings')
        settings = metadata.settings;
    else
        settings = struct('geometry',metadata.mainPara.MNZ_result, ...
            'pixelSize',metadata.mainPara.PixelSize,'wavelength',metadata.mainPara.WaveLength);
    end
    geometry = settings.geometry;
    if channel == 1
        referencePixelSize = settings.pixelSize;
    elseif settings.pixelSize ~= referencePixelSize
        error('Lensless:IncompatibleChannels', ...
            'Translation-only RGB fusion requires matching channel pixel sizes.');
    end
    wavelengths(channel) = settings.wavelength;
    channelImageCounts(channel) = numel(metadata.zPositions);
    distanceStepsByChannel{channel} = [metadata.zPositions(1), diff(metadata.zPositions(:).')];
    reconstruction = load(resultFile);
    grayscaleFolder = fullfile(outputRoot,'grayscale',channelNames{channel});
    exportGrayscaleResult(reconstruction,grayscaleFolder,grayscaleOutput,resultFile);
    channelAmplitude{channel} = double(gatherIfNeeded(abs(reconstruction.object)));
    channelMask{channel} = logical(gatherIfNeeded(reconstruction.validMask));
    channelCoverage{channel} = double(gatherIfNeeded(reconstruction.sampleCoverage));
    channelGeometry{channel} = geometry;
    focusDistances(channel) = reconstruction.focusDistance;
    channelRunFolders{channel} = runFolder;
end

%% Register wavelength channels and their masks together
% Reproduce the old display convention before any crop. Its padded-field
% extrema deliberately compress wavelength differences in the final FOV.
legacyChannels = cellfun(@mat2gray, channelAmplitude, ...
    'UniformOutput', false);
reference = 1;
referenceGeometry = channelGeometry{reference};
referenceKx = referenceGeometry.orig_M / referenceGeometry.Z;
referenceKy = referenceGeometry.orig_N / referenceGeometry.Z;
referenceDistance = focusDistances(reference);
channelShifts = zeros(channelCount, 2);

for channel = 1:channelCount
    geometry = channelGeometry{channel};
    shiftX = focusDistances(channel) * geometry.orig_M / geometry.Z - ...
        referenceDistance * referenceKx;
    shiftY = focusDistances(channel) * geometry.orig_N / geometry.Z - ...
        referenceDistance * referenceKy;
    channelShifts(channel, :) = [shiftX, shiftY];

    channelAmplitude{channel} = imtranslate(channelAmplitude{channel}, ...
        [shiftX, shiftY], 'cubic', 'OutputView', 'same', 'FillValues', 0);
    legacyChannels{channel} = imtranslate(legacyChannels{channel}, ...
        [shiftX, shiftY], 'cubic', 'OutputView', 'same', 'FillValues', 0);
    shiftedMask = imtranslate(double(channelMask{channel}), ...
        [shiftX, shiftY], 'nearest', 'OutputView', 'same', 'FillValues', 0);
    channelMask{channel} = shiftedMask >= 0.5;
    channelCoverage{channel} = imtranslate(channelCoverage{channel}, ...
        [shiftX, shiftY], 'nearest', 'OutputView', 'same', 'FillValues', 0);
end

%% Find the common physical field of view and crop every channel identically
[physicalBounds, cropBounds] = commonPhysicalBounds( ...
    channelMask, channelGeometry);
for channel = 1:channelCount
    bounds = cropBounds(channel, :);
    channelAmplitude{channel} = channelAmplitude{channel}( ...
        bounds(2):bounds(4), bounds(1):bounds(3));
    legacyChannels{channel} = legacyChannels{channel}( ...
        bounds(2):bounds(4), bounds(1):bounds(3));
    channelMask{channel} = channelMask{channel}( ...
        bounds(2):bounds(4), bounds(1):bounds(3));
    channelCoverage{channel} = channelCoverage{channel}( ...
        bounds(2):bounds(4), bounds(1):bounds(3));
end

% Refine the MNZ prediction with translation-only phase correlation. This
% corrects residual calibration and focus errors without introducing scale
% or rotation changes that would distort the reconstructed object.
residualShifts = zeros(channelCount, 2);
registrationReference = normalizeChannel( ...
    channelAmplitude{reference}, channelMask{reference});
for channel = 1:channelCount
    if channel == reference, continue; end
    registrationMoving = normalizeChannel( ...
        channelAmplitude{channel}, channelMask{channel});
    transform = imregcorr( ...
        registrationMoving, registrationReference, 'translation');
    residualShift = [transform.T(3, 1), transform.T(3, 2)];
    if any(abs(residualShift) > 25)
        warning(['Ignoring implausible residual shift for channel %s: ', ...
            '[%.3f, %.3f] pixels.'], channelNames{channel}, residualShift);
        continue;
    end
    residualShifts(channel, :) = residualShift;
    channelAmplitude{channel} = imtranslate(channelAmplitude{channel}, ...
        residualShift, 'cubic', 'OutputView', 'same', 'FillValues', 0);
    legacyChannels{channel} = imtranslate(legacyChannels{channel}, ...
        residualShift, 'cubic', 'OutputView', 'same', 'FillValues', 0);
    shiftedMask = imtranslate(double(channelMask{channel}), ...
        residualShift, 'nearest', 'OutputView', 'same', 'FillValues', 0);
    channelMask{channel} = shiftedMask >= 0.5;
    channelCoverage{channel} = imtranslate(channelCoverage{channel}, ...
        residualShift, 'nearest', 'OutputView', 'same', 'FillValues', 0);
end
channelShifts = channelShifts + residualShifts;

commonMask = all(cat(3, channelMask{:}), 3);
if ~all(commonMask(:))
    % Subpixel registration can leave a one-pixel disagreement at a border.
    [rows, columns] = find(commonMask);
    if isempty(rows)
        error('RGB registration produced an empty common valid mask.');
    end
    innerBounds = [min(columns), min(rows), max(columns), max(rows)];
    for channel = 1:channelCount
        channelAmplitude{channel} = cropByBounds( ...
            channelAmplitude{channel}, innerBounds);
        legacyChannels{channel} = cropByBounds( ...
            legacyChannels{channel}, innerBounds);
        channelMask{channel} = cropByBounds( ...
            channelMask{channel}, innerBounds);
        channelCoverage{channel} = cropByBounds( ...
            channelCoverage{channel}, innerBounds);
    end
    commonMask = all(cat(3, channelMask{:}), 3);
    physicalBounds = physicalBounds([1,2,1,2]) + innerBounds - 1;
    cropBounds = cropBounds(:,[1,2,1,2]) + innerBounds - 1;
end
commonCoverage = min(cat(3, channelCoverage{:}), [], 3);
commonCoverage(~commonMask) = 0;

%% Channel normalization and RGB composition
% Keep the former independent stretch as a diagnostic. It gives every
% wavelength its own black/white points and can therefore exaggerate small
% inter-channel differences into strong false colour.
normalizedChannels = cell(1, channelCount);
for channel = 1:channelCount
    normalizedChannels{channel} = normalizeChannel( ...
        channelAmplitude{channel}, commonMask);
    writeMaskedPNG(normalizedChannels{channel}, fullfile( ...
        outputRoot, ['amplitude_' channelNames{channel} '.png']), commonMask);
end

colorImageUncorrected = cat(3, normalizedChannels{rgbOrder(1)}, ...
    normalizedChannels{rgbOrder(2)}, normalizedChannels{rgbOrder(3)});
colorImageUncorrected(~repmat(commonMask, 1, 1, 3)) = 0;

legacyImage = cat(3, legacyChannels{rgbOrder(1)}, ...
    legacyChannels{rgbOrder(2)}, legacyChannels{rgbOrder(3)});
legacyImage(~repmat(commonMask, 1, 1, 3)) = 0;
legacyContrastImage = enhanceLuminanceContrast(legacyImage, commonMask);

% Preferred linear baseline: match only one robust brightness number per
% reconstructed wavelength, then apply the same black/white mapping to all
% three channels. This retains high contrast without forcing each channel
% independently through the complete [0, 1] range.
[linearChannels, reconstructionBrightnessGains, sharedLimits] = ...
    normalizeChannelsShared(channelAmplitude, commonMask);

if colorBalance.enabled
    [balancedChannels, balanceGainMaps] = balanceLowFrequencyShading( ...
        normalizedChannels, commonMask, colorBalance);
else
    balancedChannels = normalizedChannels;
    balanceGainMaps = ones([size(commonMask), channelCount]);
end
if colorBalance.flattenLowCoverage
    availableCoverage = min(channelImageCounts);
    fullColorCoverage = max(options.output.minimumCoverageCount + 1, ...
        ceil(colorBalance.fullColorCoverageFraction * availableCoverage));
    fullColorCoverage = min(fullColorCoverage, availableCoverage);
    colorBalance.effectiveFullColorCoverageCount = fullColorCoverage;
    [balancedChannels, lowCoverageNeutralWeight] = ...
        neutralizeLowCoverage(balancedChannels, commonMask, ...
        commonCoverage, options.output.minimumCoverageCount, ...
        fullColorCoverage);
else
    lowCoverageNeutralWeight = zeros(size(commonMask));
end
colorImageLocallyBalanced = cat(3, balancedChannels{rgbOrder(1)}, ...
    balancedChannels{rgbOrder(2)}, balancedChannels{rgbOrder(3)});
colorImageLocallyBalanced(~repmat(commonMask, 1, 1, 3)) = 0;

% Apply only the coverage-confidence correction to the shared-scale result.
% Local multiplicative shading balance is intentionally excluded here so
% that this image remains a clean test of brightness normalization.
linearDisplayChannels = linearChannels;
if colorBalance.flattenLowCoverage
    [linearDisplayChannels, linearLowCoverageNeutralWeight] = ...
        neutralizeLowCoverage(linearDisplayChannels, commonMask, ...
        commonCoverage, options.output.minimumCoverageCount, ...
        fullColorCoverage);
else
    linearLowCoverageNeutralWeight = zeros(size(commonMask));
end
colorImageLinear = cat(3, linearDisplayChannels{rgbOrder(1)}, ...
    linearDisplayChannels{rgbOrder(2)}, linearDisplayChannels{rgbOrder(3)});
colorImageLinear(~repmat(commonMask, 1, 1, 3)) = 0;

% Primary output 1: intensity-domain RGB with a global white reference and
% no coverage taper. This reproduces the former comparison image named
% 00_physical_before_coverage_taper.png.
channelIntensity = cellfun(@(amplitude) double(amplitude).^2, ...
    channelAmplitude, 'UniformOutput', false);
[physicalChannels, physicalWhiteGains, physicalWhiteLevels, ...
        physicalWhiteMask, physicalWhiteCoverageThreshold] = ...
    composePhysicalIntensity(channelIntensity, commonMask, ...
    commonCoverage, physicalIntensity);
colorImagePhysical = cat(3, physicalChannels{rgbOrder(1)}, ...
    physicalChannels{rgbOrder(2)}, physicalChannels{rgbOrder(3)});
colorImagePhysical(~repmat(commonMask, 1, 1, 3)) = 0;

% Default output: retain the old, stable colour convention and enhance only
% luminance. The more aggressive alternatives remain in comparisons/.
colorImage = legacyContrastImage;

% Smooth chroma only; preserve luminance detail from the reconstruction.
ycbcrImage = rgb2ycbcr(colorImage);
chromaKernel = fspecial('average', 2);
ycbcrImage(:, :, 2) = imfilter( ...
    ycbcrImage(:, :, 2), chromaKernel, 'replicate');
ycbcrImage(:, :, 3) = imfilter( ...
    ycbcrImage(:, :, 3), chromaKernel, 'replicate');
colorImageChromaSmoothed = ycbcr2rgb(ycbcrImage);

comparisonFolder = fullfile(outputRoot, 'comparisons');
if ~isfolder(comparisonFolder), mkdir(comparisonFolder); end
writeMaskedPNG(colorImageUncorrected, fullfile(comparisonFolder, ...
    '01_independent_percentile.png'), commonMask);
writeMaskedPNG(colorImageLocallyBalanced, fullfile(comparisonFolder, ...
    '02_independent_with_local_balance.png'), commonMask);
writeMaskedPNG(colorImageLinear, fullfile(comparisonFolder, ...
    '03_input_equalized_shared_scale.png'), commonMask);
% The two selected methods are primary outputs; all other variants remain
% under comparisons/.
writeMaskedPNG(colorImagePhysical, fullfile(outputRoot, ...
    'color_fusion_physical.png'), commonMask);
writeMaskedPNG(legacyContrastImage, fullfile(outputRoot, ...
    'color_fusion_legacy.png'), commonMask);
% Keep the legacy image at the historical filename for compatibility.
writeMaskedPNG(colorImage, fullfile(outputRoot, 'color_fusion_result.png'), commonMask);
writeMaskedPNG(colorImageChromaSmoothed, fullfile( ...
    comparisonFolder, '02_independent_with_local_balance_yuv.png'), commonMask);
writeMaskedPNG(legacyImage, fullfile(comparisonFolder, ...
    '04_legacy_full_field_scale.png'), commonMask);
save(fullfile(outputRoot, 'color_fusion_result.mat'), ...
    'colorImage', 'colorImagePhysical', 'physicalIntensity', ...
    'physicalWhiteGains', 'physicalWhiteLevels', 'physicalWhiteMask', ...
    'physicalWhiteCoverageThreshold', ...
    'colorImageLinear', 'colorImageLocallyBalanced', ...
    'colorImageUncorrected', 'colorImageChromaSmoothed', ...
    'legacyImage', 'legacyContrastImage', ...
    'commonMask', 'colorBalance', 'balanceGainMaps', ...
    'commonCoverage', 'lowCoverageNeutralWeight', ...
    'linearLowCoverageNeutralWeight', 'inputBrightness', ...
    'inputBrightnessGains', 'inputBrightnessLevels', ...
    'reconstructionBrightnessGains', 'sharedLimits', ...
    'channelShifts', 'residualShifts', 'focusDistances', 'physicalBounds', ...
    'channelRunFolders', 'channelFolders', 'channelNames', ...
    'channelImageCounts', 'wavelengths', 'distanceStepsByChannel');

if options.showFigures
    figure('Color', 'w', 'Name', 'Color APRW reconstruction');
    subplot(1, 2, 1);
    imshow(colorImagePhysical);
    title('Physical intensity with global white reference');
    subplot(1, 2, 2);
    imshow(legacyContrastImage);
    title('Legacy colour with luminance contrast');
end
fprintf('\nColor reconstruction saved to:\n  %s\n', outputRoot);

end

function [resultFile, iterationFolder] = latestReconstructionFile(runFolder)
iterations = [dir(fullfile(runFolder, 'iteration_*')); dir(fullfile(runFolder, 'Iter_*'))];
iterations = iterations([iterations.isdir]);
if isempty(iterations)
    error('No recorded APRW iteration exists in: %s', runFolder);
end
[~, order] = sort({iterations.name});
latestFolder = iterations(order(end)).name;
iterationFolder = fullfile(runFolder, latestFolder);
resultFile = fullfile( ...
    iterationFolder, 'results', 'reconstruction.mat');
if ~isfile(resultFile)
    error('Reconstruction result does not exist: %s', resultFile);
end
end

function runFolder = latestCompatibleRunFolder(channelOutputRoot, sourceInfo)
runFolder = '';
if ~isfolder(channelOutputRoot), return; end
runs = [dir(fullfile(channelOutputRoot, 'reconstruction_*')); dir(fullfile(channelOutputRoot, 'APRW_*'))];
runs = runs([runs.isdir]);
if isempty(runs), return; end
[~, order] = sort({runs.name});
for index = numel(order):-1:1
    candidate = fullfile( ...
        runs(order(index)).folder, runs(order(index)).name);
    sourceFile = fullfile(candidate, 'color_source.mat');
    if ~isfile(sourceFile), continue; end
    cached = load(sourceFile, 'sourceInfo');
    if ~isfield(cached, 'sourceInfo') || ...
            ~isequaln(cached.sourceInfo, sourceInfo)
        continue;
    end
    try
        latestReconstructionFile(candidate);
        runFolder = candidate;
        return;
    catch
        % Ignore incomplete runs and continue searching older caches.
    end
end
end

function sourceInfo = reconstructionSourceInfo(inputFile, options)
% A prepared channel MAT is the public boundary: no raw dataset paths needed.
entry = dir(inputFile);
sourceInfo.version = 6; % Full canvas plus separate shared grayscale rendering
sourceInfo.inputFile = char(inputFile);
sourceInfo.inputBytes = entry.bytes;
sourceInfo.inputModified = entry.datenum;
sourceInfo.options = options;
end

function [physicalBounds, cropBounds] = commonPhysicalBounds(masks, geometries)
channelCount = numel(masks);
minimumX = -inf;
maximumX = inf;
minimumY = -inf;
maximumY = inf;
for channel = 1:channelCount
    [rows, columns] = find(masks{channel});
    if isempty(rows)
        error('Channel %d has an empty valid-field mask.', channel);
    end
    padY = geometries{channel}.padSize(1);
    padX = geometries{channel}.padSize(2);
    minimumX = max(minimumX, min(columns) - padX);
    maximumX = min(maximumX, max(columns) - padX);
    minimumY = max(minimumY, min(rows) - padY);
    maximumY = min(maximumY, max(rows) - padY);
end
if minimumX > maximumX || minimumY > maximumY
    error('The wavelength channels have no common valid field of view.');
end

physicalBounds = [minimumX, minimumY, maximumX, maximumY];
cropBounds = zeros(channelCount, 4);
for channel = 1:channelCount
    padY = geometries{channel}.padSize(1);
    padX = geometries{channel}.padSize(2);
    cropBounds(channel, :) = [minimumX + padX, minimumY + padY, ...
        maximumX + padX, maximumY + padY];
end
end

function output = cropByBounds(input, bounds)
output = input(bounds(2):bounds(4), bounds(1):bounds(3));
end

function output = normalizeChannel(input, mask)
values = input(mask & isfinite(input));
if isempty(values)
    error('A color channel contains no finite valid pixels.');
end
limits = prctile(values, [1, 99]);
output = (input - limits(1)) / max(limits(2) - limits(1), eps);
output = min(max(output, 0), 1);
output(~mask) = 0;
end

function [images, gains, levels] = equalizeMeasurementBrightness( ...
        images, masks, settings)
% Match the robust scalar brightness of frames within one wavelength.
count = numel(images);
levels = zeros(1, count);
for index = 1:count
    values = double(images{index}(masks{index} & isfinite(images{index})));
    limits = prctile(values, settings.percentiles);
    trimmed = values(values >= limits(1) & values <= limits(2));
    levels(index) = mean(trimmed, 'omitnan');
end
positiveLevels = levels(isfinite(levels) & levels > 0);
if isempty(positiveLevels)
    error('Input frames contain no positive finite brightness values.');
end
targetLevel = median(positiveLevels);
gains = targetLevel ./ max(levels, eps);
minimumGain = 1 / settings.maximumGain;
gains = min(max(gains, minimumGain), settings.maximumGain);
for index = 1:count
    images{index} = images{index} .* gains(index);
end
end

function [channels, gains, sharedLimits] = normalizeChannelsShared( ...
        amplitudes, mask)
% Match one robust channel level, then use common display limits for RGB.
channelCount = numel(amplitudes);
levels = zeros(1, channelCount);
for channel = 1:channelCount
    values = double(amplitudes{channel}( ...
        mask & isfinite(amplitudes{channel})));
    levels(channel) = median(values, 'omitnan');
end
targetLevel = exp(mean(log(max(levels, eps))));
gains = targetLevel ./ max(levels, eps);

adjusted = cell(1, channelCount);
allValues = [];
for channel = 1:channelCount
    adjusted{channel} = double(amplitudes{channel}) .* gains(channel);
    values = adjusted{channel}(mask & isfinite(adjusted{channel}));
    allValues = [allValues; values]; %#ok<AGROW>
end
sharedLimits = prctile(allValues, [1, 99]);
denominator = max(sharedLimits(2) - sharedLimits(1), eps);
channels = cell(1, channelCount);
for channel = 1:channelCount
    output = (adjusted{channel} - sharedLimits(1)) / denominator;
    output = min(max(output, 0), 1);
    output(~mask) = 0;
    channels{channel} = output;
end
end

function output = enhanceLuminanceContrast(input, mask)
% Increase only luminance contrast; leave the legacy chroma coordinates.
ycbcr = rgb2ycbcr(input);
luminance = ycbcr(:, :, 1);
limits = prctile(luminance(mask & isfinite(luminance)), [1, 99]);
luminance = (luminance - limits(1)) / ...
    max(limits(2) - limits(1), eps);
ycbcr(:, :, 1) = min(max(luminance, 0), 1);
output = min(max(ycbcr2rgb(ycbcr), 0), 1);
output(~repmat(mask, 1, 1, 3)) = 0;
end

function [balanced, gainMaps] = balanceLowFrequencyShading( ...
        channels, mask, settings)
% Equalize broad RGB shading without changing the reconstructed channels.
channelCount = numel(channels);
maskWeight = double(mask);
sigma = settings.backgroundSigma;
normalizer = imgaussfilt(maskWeight, sigma, 'Padding', 'replicate');
background = zeros([size(mask), channelCount]);
for channel = 1:channelCount
    numerator = imgaussfilt( ...
        double(channels{channel}) .* maskWeight, sigma, ...
        'Padding', 'replicate');
    background(:, :, channel) = numerator ./ max(normalizer, eps);
end

% A geometric mean gives the three wavelengths equal influence and avoids
% selecting one possibly shaded channel as the reference.
targetBackground = exp(mean(log(max(background, eps)), 3));
minimumGain = 1 / settings.maximumGain;
maximumGain = settings.maximumGain;
gainMaps = targetBackground ./ max(background, eps);
gainMaps = min(max(gainMaps, minimumGain), maximumGain);

balanced = cell(size(channels));
for channel = 1:channelCount
    corrected = double(channels{channel}) .* gainMaps(:, :, channel);
    corrected = min(max(corrected, 0), 1);
    corrected(~mask) = 0;
    balanced{channel} = corrected;
end
end

function [channels, neutralWeight] = neutralizeLowCoverage( ...
        channels, mask, coverage, minimumCoverage, fullColorCoverage)
% Force unreliable edge chroma toward gray while retaining its luminance.
denominator = max(fullColorCoverage - minimumCoverage, 1);
colorConfidence = (coverage - minimumCoverage) / denominator;
colorConfidence = min(max(colorConfidence, 0), 1);
% Smoothstep prevents a visible boundary at the end of the transition.
colorConfidence = colorConfidence.^2 .* (3 - 2 * colorConfidence);
neutralWeight = 1 - colorConfidence;
neutralWeight(~mask) = 0;

neutral = mean(cat(3, channels{:}), 3);
for channel = 1:numel(channels)
    channels{channel} = (1 - neutralWeight) .* channels{channel} + ...
        neutralWeight .* neutral;
    channels{channel}(~mask) = 0;
end
end

function [channels, gains, whiteLevels, whiteMask, coverageThreshold] = ...
        composePhysicalIntensity(intensities, mask, coverage, settings)
% Convert three reconstructed intensity channels to display RGB using only
% one scalar white-reference gain per wavelength. The same zero point is
% retained for every channel, and no spatial or coverage taper is applied.
channelCount = numel(intensities);
maximumCoverage = max(coverage(mask));
coverageThreshold = max(1, ...
    ceil(settings.minimumCoverageFraction * maximumCoverage));
reliableMask = mask & coverage >= coverageThreshold;
if nnz(reliableMask) < settings.minimumWhitePixels
    reliableMask = mask;
end

relativeBrightness = zeros([size(mask), channelCount]);
for channel = 1:channelCount
    input = max(double(intensities{channel}), 0);
    values = input(reliableMask & isfinite(input));
    if isempty(values)
        error('Channel %d has no finite pixels for white balance.', channel);
    end
    scale = prctile(values, 99);
    relativeBrightness(:, :, channel) = input ./ max(scale, eps);
end

% A reference pixel must be bright in all three wavelengths. Percentiles
% locate the shared bright region only; they are not separate RGB stretches.
jointBrightness = min(relativeBrightness, [], 3);
threshold = prctile(jointBrightness(reliableMask), ...
    settings.whiteReferencePercentile);
whiteMask = reliableMask & jointBrightness >= threshold;
if nnz(whiteMask) < settings.minimumWhitePixels
    retainedFraction = min(1, ...
        settings.minimumWhitePixels / nnz(reliableMask));
    threshold = prctile(jointBrightness(reliableMask), ...
        100 * (1 - retainedFraction));
    whiteMask = reliableMask & jointBrightness >= threshold;
end

whiteLevels = zeros(1, channelCount);
gains = zeros(1, channelCount);
channels = cell(1, channelCount);
for channel = 1:channelCount
    input = max(double(intensities{channel}), 0);
    whiteLevels(channel) = median(input(whiteMask), 'omitnan');
    gains(channel) = settings.targetWhite / max(whiteLevels(channel), eps);
    output = min(max(input .* gains(channel), 0), 1);
    output(~mask) = 0;
    channels{channel} = output;
end
end
