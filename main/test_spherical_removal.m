function outputFolder = test_spherical_removal(iterationFolder)
%TEST_SPHERICAL_REMOVAL Compare identical saved fields before TV/cropping.
% Call with no argument to use the most recently saved reconstruction, or
% pass a specific Iter_XXXX folder. No APRW iterations are repeated.
projectRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(projectRoot, 'src'));
if nargin < 1 || isempty(iterationFolder)
    files = dir(fullfile(projectRoot, 'ResultFolder', 'APRW_*', ...
        'Iter_*', 'results', 'reconstruction.mat'));
    if isempty(files), error('No saved reconstruction found.'); end
    [~, index] = max([files.datenum]);
    iterationFolder = fileparts(files(index).folder);
end
reconstruction = load(fullfile(iterationFolder, 'results', ...
    'reconstruction.mat'), 'field', 'focusDistance', 'bounds');
metadata = load(fullfile(iterationFolder, 'diagnostics', 'meta.mat'), 'mainPara');
parameters = metadata.mainPara;
referenceField = reconstruction.field;
referenceIllumination = parameters.IllumSet;
assert(isequal(size(referenceField), size(referenceIllumination)), ...
    'Saved illumination and field must have identical dimensions.');
focusDistance = reconstruction.focusDistance;

% Remove only the saved reference-plane illumination phase. Do not divide
% by its amplitude or remove the sample-plane sphere a second time.
flattenedField = referenceField .* conj(referenceIllumination ./ ...
    max(abs(referenceIllumination), eps));
baselineField = propagate(referenceField, parameters.PixelSize, ...
    parameters.WaveLength, -focusDistance);
testField = propagate(flattenedField, parameters.PixelSize, ...
    parameters.WaveLength, -focusDistance);
baselineAmplitude = abs(baselineField);
testAmplitude = abs(testField);

% Same distance, pixel pitch, canvas and display limits in both branches.
% Baseline amplitude equals the normal pipeline before TV: removing the
% sample-plane phase afterwards cannot change its amplitude.
bounds = reconstruction.bounds;
rows = bounds(2):bounds(4);
columns = bounds(1):bounds(3);
baselineCrop = baselineAmplitude(rows, columns);
displayLimits = prctile(double(baselineCrop(:)), [1, 99]);
scale = max(diff(displayLimits), eps);
baselinePNG = min(max((double(baselineAmplitude)-displayLimits(1))/scale, 0), 1);
testPNG = min(max((double(testAmplitude)-displayLimits(1))/scale, 0), 1);
outputFolder = fullfile(iterationFolder, 'sphere_removal_test', ...
    char(datetime('now', 'Format', 'yyyyMMdd_HHmmss_SSS')));
mkdir(outputFolder);
imwrite(baselinePNG, fullfile(outputFolder, 'baseline_full.png'));
imwrite(testPNG, fullfile(outputFolder, 'removed_before_backprop_full.png'));
imwrite(baselinePNG(rows, columns), ...
    fullfile(outputFolder, 'baseline_same_crop.png'));
imwrite(testPNG(rows, columns), ...
    fullfile(outputFolder, 'removed_before_backprop_same_crop.png'));
save(fullfile(outputFolder, 'comparison.mat'), 'baselineField', ...
    'testField', 'focusDistance', 'displayLimits', 'bounds', 'iterationFolder', '-v7.3');
fprintf('Baseline and phase-removal images are saved separately.\n');
fprintf('Fixed distance %.6g m; no TV, amplitude cap, or support masking.\n', focusDistance);
fprintf('Comparison saved to: %s\n', outputFolder);
end
