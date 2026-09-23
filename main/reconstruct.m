clc;
close all;
clear;

% Standard reconstruction entry point.
projectRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(fullfile(projectRoot, 'src')));

inputFile = fullfile(projectRoot, 'data', 'reconstruction_input.mat');
if ~isfile(inputFile)
    error('Input data file is missing. Run create_input_data first:\n%s', inputFile);
end
loaded = load(inputFile, 'inputData');

options.iterations = 18;
options.recordEvery = 1;
options.focus.prior = 1.68e-3;
options.focus.halfRange = 3e-3;
options.focus.step = 0.01e-3;
options.tv.enabled = true;
options.tv.lambdaMin = 2e-3;  % maximum-overlap protected region
options.tv.lambdaMax = 2e-2;  % single-measurement region
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
options.output.validMaskThreshold = 0.01;
% Default PNG: largest offset-expanded FOV with multi-plane support.
% The exact first-camera FOV remains available as originalFOV_amplitude.png.
options.output.defaultToOriginalFOV = false;
% Retain the complete offset-expanded union. Sample-plane displacement is
% handled later when the reconstruction is cropped after back propagation.
options.output.minimumCoverageCount = 1;
options.output.limitHorizontalToReferenceFOV = false;
options.output.trimPropagationBoundary = false;
options.showFigures = true;

result = APRW(loaded.inputData, options);
