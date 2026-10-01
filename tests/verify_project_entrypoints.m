function verify_project_entrypoints()
%VERIFY_PROJECT_ENTRYPOINTS Exercise public workflows without private adapters.
projectRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(projectRoot);
setup_project;
testRoot = fullfile(projectRoot, 'outputs', 'entrypoint_verification', ...
    char(datetime('now', 'Format', 'yyyyMMdd_HHmmss_SSS')));
mkdir(testRoot);
inputData = syntheticInput();
inputFile = fullfile(testRoot, 'input.mat');
save(inputFile, 'inputData');
initialFigures = findall(groot, 'Type', 'figure');
for name = {'reconstruct', 'reconstruct_fast'}
    cfg = loadProjectConfig(name{1}, false);
    cfg.inputFile = inputFile;
    cfg.options = smallOptions(cfg.options);
    cfg.options.focus.halfRange = 0.04e-3;
    cfg.options.focus.step = 0.02e-3;
    if strcmp(name{1},'reconstruct_fast'), cfg.options.focus.method = 'adfrft'; end
    cfg.options.output.rootFolder = fullfile(testRoot, name{1});
    [fields, runFolder] = feval(name{1}, cfg);
    assert(~isempty(fields));
    iterationFolder = fullfile(runFolder, 'iteration_0001');
    focus = load(fullfile(iterationFolder, 'diagnostics', 'autofocus.mat'));
    result = load(fullfile(iterationFolder, 'results', 'reconstruction.mat'));
    assert(numel(focus.distances) == 5 && focus.searchPerformed);
    assert(focus.bestDistance == focus.distances(focus.bestIndex));
    assert(focus.focusMetric(focus.bestIndex) == max(focus.focusMetric));
    assert(result.focusDistance == focus.bestDistance);
    assert(strcmp(focus.method,cfg.options.focus.method));
    assert(isfile(fullfile(iterationFolder, 'diagnostics', 'autofocus.png')));
    assert(isequal(findall(groot, 'Type', 'figure'), initialFigures));
    % The plotted values are actual saved scores with metres converted to mm.
    amplitude = imread(fullfile(iterationFolder, 'results', 'amplitude.png'));
    phase = imread(fullfile(iterationFolder, 'results', 'phase_heatmap.png'));
    fig = plotAutofocus(focus, amplitude, phase, 'Regression', [], 'off');
    cleanup = onCleanup(@() close(fig));
    curve = findobj(fig, 'Tag', 'FocusCurve');
    assert(isequal(curve.XData(:), focus.distances(:)*1e3));
    assert(isequal(curve.YData(:), focus.focusMetric(:)));
    clear cleanup;
end

cfg = loadProjectConfig('color', false);
cfg.inputFile = fullfile(testRoot, 'color_input.mat');
cfg.outputRoot = fullfile(testRoot, 'color_output');
cfg.options = smallOptions(cfg.options);
cfg.reuseExistingReconstructions = true;
colorInput.channelFiles = cell(1,3);
wavelengths = [460e-9, 514e-9, 647e-9];
for channel = 1:3
    inputData = syntheticInput(wavelengths(channel));
    colorInput.channelFiles{channel} = sprintf('channel_%d.mat', channel);
    save(fullfile(testRoot, colorInput.channelFiles{channel}), 'inputData');
end
save(cfg.inputFile, 'colorInput');
outputRoot = color(cfg);
first = load(fullfile(outputRoot, 'color_fusion_result.mat'));
assert(size(first.colorImage,3)==3 && all(isfinite(first.colorImage(:))));
assert(any(first.commonMask(:)));
color(cfg);
second = load(fullfile(outputRoot, 'color_fusion_result.mat'));
assert(isequal(first.channelRunFolders, second.channelRunFolders));
assert(isequal(first.colorImage, second.colorImage));
assert(isequal(findall(groot, 'Type', 'figure'), initialFigures));
fprintf('Public entry point, autofocus plot and color cache checks passed.\n');
end

function options = smallOptions(options)
options.iterations = 1;
options.recordEvery = 1;
options.focus.prior = 1e-3;
options.focus.halfRange = 0;
options.tv.enabled = false;
options.showFigures = false;
end
