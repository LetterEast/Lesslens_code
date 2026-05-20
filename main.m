clc,close all,clear
addpath(genpath(pwd));
%==========================================================================
% Multi-distance lensless reconstruction pipeline (APRW)
%
% Steps:
%   1) Raw image averaging: RawImgAverage() -> RawImgSet
%   2) Load averaged images: loadImg() -> img_set (cell, multi-plane intensity)
%   3) Build axial interval set: DistanceIntervalSet
%   4) Calibrate spherical illumination: XYZSphereIllumCalibration() -> MNZ_result
%   5) Generate illumination field: getMultiAngleIllum() -> IllumPaSet
%   6) Iterative reconstruction + recording: getRec_APRW() -> Rec_nInterative
%
% Inputs:
%   RawImgFolder (png frames), WaveLength [m], PixelSize [m/pixel],
%   DistanceInterval [m], numImages
%
% Output:
%   Results saved by getRec_APRW to .\ResultFolder
%
% Unit convention:
%   wavelength [m], distances [m], pixel size [m/pixel]
%==========================================================================

%% Parameters setting (physical + acquisition)
RawImgFolder = "D:\Desktop\data\2026\rgb\img\02_green_sample";
load('D:\Desktop\data\2026\rgb\MNZ\02\MNZ_result.mat')

WaveLength          = 514e-9;                          % [m] wavelength
PixelSize           = 3e-6;                            % [m/pixel] sensor effective pixel size
nIterative          = 10;                            % iteration count
numImages           = 5;                               % number of captured positions/images
DistanceIntervalSet = [0,ones(1,numImages-1)]*1e-3;  % [m] distance interval between adjacent captures

%% Distance interval configuration (calibration search ranges)
% Options
iIte_record = 10;                             % record every iIte_record iterations

% CCD axial translation per capture
ParaDisInterval.DisIntervalPre =  0.5e-3;      % [m] prior interval
ParaDisInterval.DisIntervalHalfRange = 0.5e-3;% [m] search half range
ParaDisInterval.rough = 0.01e-3;            % [m] coarse step

%  Sample-to-CCD distance
ParaD_Sample2CCD.D_Sample2CCDPre       = 2e-3;     % [m]
ParaD_Sample2CCD.D_Sample2CCDHalfRange = 1e-3;     % [m]
ParaD_Sample2CCD.rough                 = 0.01e-3;  % [m] 粗测步长
ParaD_Sample2CCD.Interval              = 0.1e-3;   % [m]

MainPara.WaveLength       = WaveLength;
MainPara.PixelSize        = PixelSize;
MainPara.nIterative       = nIterative;
MainPara.ParaD_Sample2CCD = ParaD_Sample2CCD;
MainPara.MNZ_result       = MNZ_result;
%% 1) Raw image average -> RawImgSet (3D array)
OutputPath = RawImgAverage(RawImgFolder, numImages);

%% 2) Load cropped images as cell stack for reconstruction
img_set = loadImg(OutputPath);
%% 3) Distance interval set
% IllumPaSet = getMultiAngleIllum(img_set{1},MNZ_result, PixelSize, WaveLength);
% DistanceIntervalSet = DisIntervalCalibration(img_set,numImages,ParaDisInterval,PixelSize,WaveLength,IllumPaSet);
% DistanceIntervalSet = AutoFocusSphericalWave(img_set,numImages,MainPara);
%% 4) Calibrate Illumimation
[out_img_set,MNZ_result] = TiltIllumination(img_set,MNZ_result,PixelSize,DistanceIntervalSet);
IllumPaSet               = getMultiAngleIllum(out_img_set{1},MNZ_result, PixelSize, WaveLength);

%% 5) Reconstruction
MainPara.IllumSet         = IllumPaSet;  % placeholder illumination (explain why)
MainPara.MNZ_result       = MNZ_result;
Rec_nInterative = getRec_APRW(out_img_set, MainPara, DistanceIntervalSet, iIte_record);