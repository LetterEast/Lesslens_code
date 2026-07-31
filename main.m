clc;
close all;
clear;

projectRoot = fileparts(mfilename('fullpath'));
addpath(genpath(projectRoot));

%% User parameters
cfg.data.rawFolder = 'D:\Desktop\data\2026\4.20\1.0';
cfg.data.calibrationFile = 'D:\Desktop\data\2026\4.20\MNZ_result.mat';
cfg.data.numImages = 5;

cfg.optics.wavelength = 514e-9;       % [m]
cfg.optics.pixelSize  = 3e-6;         % [m/pixel]
cfg.distance.steps = [0, ones(1, cfg.data.numImages - 1)] * 1e-3;

cfg.reconstruction.iterations  = 9;
cfg.reconstruction.recordEvery = 3;
cfg.output.cropToValidFOV     = true;
cfg.output.validMaskThreshold = 0.01;
cfg.output.zeroFillInvalid    = true;

cfg.calibration.sampleToCCD.prior     = 2e-3;    % [m]
cfg.calibration.sampleToCCD.halfRange = 1.5e-3;  % [m]
cfg.calibration.sampleToCCD.step      = 0.01e-3; % [m]
cfg.calibration.sampleToCCD.interval  = 0.1e-3;  % [m]

%% Load and validate
calibration = load(cfg.data.calibrationFile, 'MNZ_result');
MNZ_result = calibration.MNZ_result;
assert(numel(cfg.distance.steps) == cfg.data.numImages, ...
    'cfg.distance.steps must contain one value per image.');

%% Pre-process and align tilted illumination
averagedFolder = RawImgAverage(cfg.data.rawFolder, cfg.data.numImages);
imgSet = loadImg(averagedFolder);
assert(numel(imgSet) == cfg.data.numImages, ...
    'Expected %d images, but loaded %d.', cfg.data.numImages, numel(imgSet));
[imgSet, MNZ_result] = TiltIllumination( ...
    imgSet, MNZ_result, cfg.optics.pixelSize, cfg.distance.steps);

%% Reconstruction parameters
sampleToCCD.D_Sample2CCDPre       = cfg.calibration.sampleToCCD.prior;
sampleToCCD.D_Sample2CCDHalfRange = cfg.calibration.sampleToCCD.halfRange;
sampleToCCD.rough                 = cfg.calibration.sampleToCCD.step;
sampleToCCD.Interval              = cfg.calibration.sampleToCCD.interval;

mainPara.WaveLength       = cfg.optics.wavelength;
mainPara.PixelSize        = cfg.optics.pixelSize;
mainPara.nIterative       = cfg.reconstruction.iterations;
mainPara.ParaD_Sample2CCD = sampleToCCD;
mainPara.MNZ_result       = MNZ_result;
mainPara.Output           = cfg.output;
mainPara.IllumSet = getMultiAngleIllum( ...
    imgSet{1}, MNZ_result, cfg.optics.pixelSize, cfg.optics.wavelength);

Rec_nInterative = getRec_APRW( ...
    imgSet, mainPara, cfg.distance.steps, cfg.reconstruction.recordEvery);
