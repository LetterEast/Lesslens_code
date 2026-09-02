function create_input_data()
%CREATE_INPUT_DATA Create the MAT file consumed by reconstruction scripts.
% Input images must already be preprocessed intensity measurements.

projectRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(fullfile(projectRoot, 'src')));

%% Dataset and acquisition geometry
imageFolder = ...
    '\\192.168.2.166\d\lesslens\2026.9.2\USAF_1951_3\G\foreground_img';
calibrationFile = fullfile(projectRoot, 'data', 'calibration', ...
    '9.1_G', 'MNZ_result.mat');
outputFile = fullfile(projectRoot, 'data', 'reconstruction_input.mat');
expectedImageCount = 61;
wavelength = 514e-9; % [m]
pixelSize = 3e-6;    % [m/pixel]
distanceSteps = [0, ones(1, expectedImageCount - 1)] * 0.1e-3; % [m]

%% Build the reconstruction input
calibration = load(calibrationFile, 'MNZ_result');
if ~isfield(calibration, 'MNZ_result')
    error('Calibration file does not contain MNZ_result: %s', calibrationFile);
end
[images, geometry] = prepareMeasurements( ...
    imageFolder, calibration.MNZ_result, distanceSteps);
if numel(images) ~= expectedImageCount
    error('Expected %d images, but found %d.', ...
        expectedImageCount, numel(images));
end

inputData = struct('images', {images}, 'geometry', geometry, ...
    'distanceSteps', distanceSteps, 'wavelength', wavelength, ...
    'pixelSize', pixelSize, 'sourceImageFolder', imageFolder, ...
    'calibrationFile', calibrationFile, 'createdAt', datetime('now'));
outputFolder = fileparts(outputFile);
if ~isfolder(outputFolder), mkdir(outputFolder); end
save(outputFile, 'inputData', '-v7.3');
fprintf('Reconstruction input saved to:\n  %s\n', outputFile);
end
