function [Object] = back2Object(ImgRec, D_Sample2CCD, MainPara)

%% Parameters setting
WaveLength          = MainPara.WaveLength;
PixelSize           = MainPara.PixelSize;

%% back to Object

validSize = size(ImgRec);
if isfield(MainPara, 'MNZ_result') && ...
        isfield(MainPara.MNZ_result, 'orig_size')
    validSize = MainPara.MNZ_result.orig_size;
end
Object = propGPU(ImgRec, PixelSize, WaveLength, ...
    -D_Sample2CCD, 0, 0, validSize); % light field propagating back on sample plane
disp('finish back to object');