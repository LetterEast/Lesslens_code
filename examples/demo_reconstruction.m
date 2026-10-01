function runFolder = demo_reconstruction(showFigures, focusMethod)
%DEMO_RECONSTRUCTION Run a small public example without personal data/tools.
if nargin < 1, showFigures = true; end
if nargin < 2, focusMethod = 'fft'; end
projectRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(projectRoot);
setup_project;
demoRoot = fullfile(projectRoot, 'outputs', 'demo', ...
    char(datetime('now','Format','yyyyMMdd_HHmmss_SSS')));
mkdir(demoRoot);
inputData = syntheticInput();
config = loadProjectConfig('reconstruct', false);
config.inputFile = fullfile(demoRoot, 'input.mat');
save(config.inputFile, 'inputData');
config.options.iterations = 3;
config.options.recordEvery = 3;
config.options.focus.prior = 1e-3;
config.options.focus.halfRange = 0.1e-3;
config.options.focus.step = 0.02e-3;
config.options.focus.method = focusMethod;
config.options.showFigures = showFigures;
config.options.output.rootFolder = demoRoot;
[~, runFolder] = reconstruct(config);
end
