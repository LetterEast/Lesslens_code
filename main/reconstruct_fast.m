clc;
close all;
clear;

% Headless reconstruction entry point.
projectRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(fullfile(projectRoot, 'src')));

inputFile = fullfile(projectRoot, 'data', 'reconstruction_input.mat');
if ~isfile(inputFile)
    error('Input data file is missing. Run create_input_data first:\n%s', inputFile);
end
loaded = load(inputFile, 'inputData');

options.iterations = 18;
options.recordEvery = options.iterations;
options.focus.prior = 1.68e-3;
options.focus.halfRange = 0e-3;
options.focus.step = 0.01e-3;
options.tv.enabled = true;
options.tv.lambdaMin = 8e-3;  % maximum-overlap protected region
options.tv.lambdaMax = 10e-2;  % single-measurement region
options.tv.coveragePower = 2;
options.tv.coverageWeight = 0.8;
options.tv.uncertaintyWeight = 0.2;
options.tv.boundaryWeight = 0.5;
options.tv.boundaryPower = 2;
options.tv.directionGain = 0.4;
options.tv.step = 2;
options.tv.subiterations = 10;
options.output.rootFolder = fullfile(projectRoot, 'ResultFolder');
options.output.cropToValidFOV = true;
options.output.zeroFillInvalid = true;
options.showFigures = false;
tic;
result = APRW(loaded.inputData, options);
toc;
