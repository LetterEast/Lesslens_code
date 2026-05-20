function [Object] = back2Object(ImgRec, D_Sample2CCD, MainPara)

%% Parameters setting
WaveLength          = MainPara.WaveLength;
PixelSize           = MainPara.PixelSize;

%% back to Object

Object = propGPU(ImgRec,PixelSize,WaveLength,-D_Sample2CCD); % light field propagating back on sample plane
disp('finish back to object');