clc;
close all;
clear;

% Color reconstruction entry point. Each wavelength uses the same APRW
% pipeline as grayscale reconstruction, followed by mask-aware RGB fusion.
projectRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(fullfile(projectRoot, 'src')));

%% Dataset configuration
imageRoot = 'D:\Desktop\data\2026\9.1\Breast_sec\UD_1';
calibrationRoot = 'D:\Desktop\data\2026\9.1\Pumpkin_sec\MNZ';
outputRoot = fullfile(projectRoot, 'ResultFolder', 'Color');
% Set COLOR_FORCE_RECONSTRUCT=1 before launching MATLAB to ignore the cache.
reuseExistingReconstructions = ...
    ~strcmp(getenv('COLOR_FORCE_RECONSTRUCT'), '1');
% Explicit order: folders 01, 02 and 03 contain B, G and R respectively.
channelFolders = {'01', '02', '03'};
channelNames = {'B', 'G', 'R'};
wavelengths = [460e-9, 514e-9, 647e-9];
rgbOrder = [3, 2, 1];

% Correct only broad channel-dependent shading during RGB fusion. The
% reconstructed grayscale channels themselves are left unchanged.
colorBalance.enabled = true;
colorBalance.backgroundSigma = 100;
colorBalance.maximumGain = 1.25;
colorBalance.flattenLowCoverage = true;
colorBalance.fullColorCoverageFraction = 0.60;

% Second primary output: convert reconstructed field amplitude to
% transmitted intensity and apply one global white gain per wavelength.
% No coverage-dependent colour taper or local spatial colour correction is
% applied to this version.
physicalIntensity.whiteReferencePercentile = 85;
physicalIntensity.minimumWhitePixels = 2048;
physicalIntensity.minimumCoverageFraction = 0.60;
physicalIntensity.targetWhite = 0.98;

% Correct scalar exposure drift between measurements before APRW. This is
% not image averaging: every frame is retained, and only one gain is
% applied to each frame. The gain limit prevents a damaged frame from
% dominating the phase-retrieval constraint.
inputBrightness.enabled = true;
inputBrightness.percentiles = [5, 95];
inputBrightness.maximumGain = 3;

pixelSize = 3e-6;
% distanceInterval = 0.1e-3;
 distanceSteps = (0:0.2:1.4) * 1e-3; % [m]
%% Shared reconstruction settings
options.iterations = 6;
options.recordEvery = options.iterations;
options.focus.prior = 1.62e-3;
options.focus.halfRange = 0.5e-3;
options.focus.step = 0.01e-3;
options.tv.enabled = true;
options.tv.lambdaMin = 2e-3;
options.tv.lambdaMax = 2e-2;
options.tv.coveragePower = 2;
options.tv.coverageWeight = 0.8;
options.tv.uncertaintyWeight = 0.2;
options.tv.boundaryWeight = 0.5;
options.tv.boundaryPower = 2;
options.tv.directionGain = 0.4;
options.tv.step = 2;
options.tv.subiterations = 10;
options.output.cropToValidFOV = false;
options.output.zeroFillInvalid = false;
options.output.validMaskThreshold = 0.01;
% Preserve each wavelength's complete offset-expanded valid FOV. RGB
% fusion below uses the intersection of the three registered masks.
options.output.minimumCoverageCount = 1;
options.output.limitHorizontalToReferenceFOV = false;
options.output.defaultToOriginalFOV = false;
options.output.trimPropagationBoundary = false;
options.showFigures = false;
options.inputBrightness = inputBrightness;

if ~isfolder(imageRoot), error('Image root does not exist: %s', imageRoot); end
if ~isfolder(calibrationRoot)
    error('Calibration root does not exist: %s', calibrationRoot);
end
if ~isfolder(outputRoot), mkdir(outputRoot); end

channelCount = numel(channelFolders);
channelAmplitude = cell(1, channelCount);
channelMask = cell(1, channelCount);
channelCoverage = cell(1, channelCount);
channelGeometry = cell(1, channelCount);
focusDistances = zeros(1, channelCount);
channelRunFolders = cell(1, channelCount);
channelImageCounts = zeros(1, channelCount);
distanceStepsByChannel = cell(1, channelCount);
inputBrightnessGains = cell(1, channelCount);
inputBrightnessLevels = cell(1, channelCount);

%% Reconstruct every wavelength with the maintained APRW implementation
for channel = 1:channelCount
    frameBrightnessGains = [];
    frameBrightnessLevels = [];
    imageFolder = fullfile(imageRoot, channelFolders{channel});
    calibrationFile = fullfile( ...
        calibrationRoot, channelFolders{channel}, 'MNZ_result.mat');
    if ~isfolder(imageFolder)
        error('Channel image folder does not exist: %s', imageFolder);
    end
    if ~isfile(calibrationFile)
        error('Channel calibration does not exist: %s', calibrationFile);
    end

    imageEntries = supportedImageEntries(imageFolder);
    channelImageCount = numel(imageEntries);
    if channelImageCount < 2
        error('Channel %s needs at least two images; found %d.', ...
            channelNames{channel}, channelImageCount);
    end
    channelImageCounts(channel) = channelImageCount;
    % distanceSteps = [0, ones(1, channelImageCount - 1)] * distanceInterval;
    
    distanceStepsByChannel{channel} = distanceSteps;

    fprintf('\n[%d/%d] Reconstructing %s channel (%g nm)\n', ...
        channel, channelCount, channelNames{channel}, ...
        wavelengths(channel) * 1e9);
    channelOutputRoot = fullfile( ...
        outputRoot, ['APRW_' channelNames{channel}]);
    options.output.rootFolder = channelOutputRoot;
    sourceInfo = reconstructionSourceInfo(imageFolder, calibrationFile, ...
        wavelengths(channel), pixelSize, distanceSteps, options);
    if reuseExistingReconstructions
        runFolder = latestCompatibleRunFolder(channelOutputRoot, sourceInfo);
    else
        runFolder = '';
    end
    if isempty(runFolder)
        calibration = load(calibrationFile, 'MNZ_result');
        if ~isfield(calibration, 'MNZ_result')
            error('MNZ_result is missing from: %s', calibrationFile);
        end
        [images, geometry] = prepareMeasurements( ...
            imageFolder, calibration.MNZ_result, distanceSteps);
        if inputBrightness.enabled
            [images, frameBrightnessGains, frameBrightnessLevels] = ...
                equalizeMeasurementBrightness(images, ...
                geometry.ValidMaskHard, inputBrightness);
            fprintf('  Input brightness gains: %.3f to %.3f\n', ...
                min(frameBrightnessGains), max(frameBrightnessGains));
        end
        inputData = struct('images', {images}, 'geometry', geometry, ...
            'distanceSteps', distanceSteps, ...
            'wavelength', wavelengths(channel), 'pixelSize', pixelSize);
        [~, runFolder] = APRW(inputData, options);
        save(fullfile(runFolder, 'color_source.mat'), 'sourceInfo', ...
            'frameBrightnessGains', 'frameBrightnessLevels');
    else
        fprintf('Reusing matching reconstruction: %s\n', runFolder);
        sourceData = load(fullfile(runFolder, 'color_source.mat'), ...
            'frameBrightnessGains', 'frameBrightnessLevels');
        if isfield(sourceData, 'frameBrightnessGains')
            frameBrightnessGains = sourceData.frameBrightnessGains;
            frameBrightnessLevels = sourceData.frameBrightnessLevels;
        end
    end
    inputBrightnessGains{channel} = frameBrightnessGains;
    inputBrightnessLevels{channel} = frameBrightnessLevels;

    [resultFile, iterationFolder] = latestReconstructionFile(runFolder);
    if ~exist('geometry', 'var')
        metadata = load(fullfile( ...
            iterationFolder, 'diagnostics', 'meta.mat'), 'mainPara');
        geometry = metadata.mainPara.MNZ_result;
    end
    reconstruction = load(resultFile, ...
        'object', 'validMask', 'focusDistance');
    channelAmplitude{channel} = double(gather(abs(reconstruction.object)));
    channelMask{channel} = logical(gather(reconstruction.validMask));
    if isfield(geometry, 'CoverageCount')
        channelCoverage{channel} = double(gather(geometry.CoverageCount));
    else
        channelCoverage{channel} = ...
            double(channelMask{channel}) * channelImageCount;
    end
    channelGeometry{channel} = geometry;
    focusDistances(channel) = reconstruction.focusDistance;
    channelRunFolders{channel} = runFolder;
    clear geometry images inputData calibration sourceInfo sourceData;
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
    imwrite(normalizedChannels{channel}, fullfile( ...
        outputRoot, ['amplitude_' channelNames{channel} '.png']));
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
imwrite(colorImageUncorrected, fullfile(comparisonFolder, ...
    '01_independent_percentile.png'));
imwrite(colorImageLocallyBalanced, fullfile(comparisonFolder, ...
    '02_independent_with_local_balance.png'));
imwrite(colorImageLinear, fullfile(comparisonFolder, ...
    '03_input_equalized_shared_scale.png'));
% The two selected methods are primary outputs; all other variants remain
% under comparisons/.
imwrite(colorImagePhysical, fullfile(outputRoot, ...
    'color_fusion_physical.png'));
imwrite(legacyContrastImage, fullfile(outputRoot, ...
    'color_fusion_legacy.png'));
% Keep the legacy image at the historical filename for compatibility.
imwrite(colorImage, fullfile(outputRoot, 'color_fusion_result.png'));
imwrite(colorImageChromaSmoothed, fullfile( ...
    comparisonFolder, '02_independent_with_local_balance_yuv.png'));
imwrite(legacyImage, fullfile(comparisonFolder, ...
    '04_legacy_full_field_scale.png'));
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

function [resultFile, iterationFolder] = latestReconstructionFile(runFolder)
iterations = dir(fullfile(runFolder, 'Iter_*'));
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
runs = dir(fullfile(channelOutputRoot, 'APRW_*'));
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

function sourceInfo = reconstructionSourceInfo( ...
        imageFolder, calibrationFile, wavelength, pixelSize, ...
        distanceSteps, options)
% Tie a cached reconstruction to its images, calibration and parameters.
entries = supportedImageEntries(imageFolder);

calibrationEntry = dir(calibrationFile);
% Increment when reconstruction/FOV semantics change so stale channel
% reconstructions are not silently reused by the color cache.
sourceInfo.version = 2;
sourceInfo.imageFolder = normalizePath(imageFolder);
sourceInfo.imageNames = lower(string({entries.name}));
sourceInfo.imageBytes = [entries.bytes];
sourceInfo.imageModified = [entries.datenum];
sourceInfo.calibrationFile = normalizePath(calibrationFile);
sourceInfo.calibrationBytes = calibrationEntry.bytes;
sourceInfo.calibrationModified = calibrationEntry.datenum;
sourceInfo.wavelength = wavelength;
sourceInfo.pixelSize = pixelSize;
sourceInfo.distanceSteps = distanceSteps;
sourceInfo.options = options;
end

function entries = supportedImageEntries(folder)
entries = dir(folder);
entries = entries(~[entries.isdir]);
validExtensions = {'.png', '.jpg', '.jpeg', '.bmp', '.tif', '.tiff'};
keep = false(size(entries));
for index = 1:numel(entries)
    [~, ~, extension] = fileparts(entries(index).name);
    keep(index) = any(strcmpi(extension, validExtensions));
end
entries = entries(keep);
[~, order] = sort(lower(string({entries.name})));
entries = entries(order);
end

function path = normalizePath(path)
path = lower(strrep(char(path), '/', '\'));
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
