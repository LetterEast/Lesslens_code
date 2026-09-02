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
options.focus.prior = 1.5e-3;
options.focus.halfRange = 0.5e-3;
options.focus.step = 0.01e-3;
options.output.rootFolder = fullfile(projectRoot, 'ResultFolder');
options.output.cropToValidFOV = true;
options.output.zeroFillInvalid = true;
options.showFigures = false;
tic;
result = APRW(loaded.inputData, options);
toc;