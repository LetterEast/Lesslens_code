function [out_img_set,MNZ_result] = TiltIllumination(img_set,MNZ_result,PixelSize,DistanceIntervalSet)

M = MNZ_result.M;
N = MNZ_result.N;
Z = MNZ_result.Z;
kx = M/Z; %X方向正切值
ky = N/Z; %Y方向正切值
title_folder = "title_folder";
if ~exist(title_folder,"dir"),mkdir(title_folder); end

num_img = length(img_set);
[mRow, nCol] = size(img_set{1});


fx = (-nCol/2 : nCol/2-1) / (nCol*PixelSize);   % cycles/m
fy = (-mRow/2 : mRow/2-1) / (mRow*PixelSize);
[u, v] = meshgrid(fx, fy);

out_img_set = cell(1,num_img);
out_img_set{1} = img_set{1};
imwrite(mat2gray(img_set{1}),fullfile(title_folder,"001.png"))
DistanceInterval = 0;
for i = 2:num_img
    Interval = DistanceIntervalSet(i);
    DistanceInterval = Interval + DistanceInterval;
    InputImg_fft = fftshift(fft2(img_set{i}));
    ShiftX = DistanceInterval*kx*PixelSize;
    ShiftY = DistanceInterval*ky*PixelSize;
    Phase = exp(-1i * 2*pi * ( u * ShiftX + v * ShiftY));
    out_img_set{i} = ifft2(ifftshift(InputImg_fft .* Phase));
    imwrite(mat2gray(abs(out_img_set{i})),fullfile(title_folder,sprintf("%0.3d.png",i)))
end

MNZ_result.M = 0;
MNZ_result.N = 0;
