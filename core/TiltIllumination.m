function [out_img_set,MNZ_result] = TiltIllumination(img_set,MNZ_result,PixelSize,DistanceIntervalSet)

M = MNZ_result.M;
N = MNZ_result.N;
Z = MNZ_result.Z;
kx = M/Z; %X方向正切值
ky = N/Z; %Y方向正切值
title_folder = "tilt_folder";
if ~exist(title_folder,"dir"),mkdir(title_folder); end

num_img = length(img_set);
[mRow_orig, nCol_orig] = size(img_set{1});

%% 固定小填充（20 像素）用于 propGPU 边界平滑
% 不再需要大的平移补偿填充，CovMask 已经负责识别有效视场区域
padSize = [20, 20];

fprintf('[TiltIllumination] padSize: [%d, %d] (fixed boundary padding)\n', padSize(1), padSize(2));

%% 对所有输入图像进行相同的边界复制填充
for i = 1:num_img
    img_set{i} = padarray(img_set{i}, padSize, 'replicate');
end

[mRow, nCol] = size(img_set{1});

%% 产生平滑边缘蒙版 W
% 在原图区域为 1，在填充区域平滑衰减到 0，消除移入的卷边和边界阶跃
maskY = ones(mRow, 1);
maskY(1:padSize(1)) = sin((0:padSize(1)-1)'/padSize(1) * pi/2).^2;
maskY(end-padSize(1)+1:end) = cos((0:padSize(1)-1)'/padSize(1) * pi/2).^2;

maskX = ones(1, nCol);
maskX(1:padSize(2)) = sin((0:padSize(2)-1)/padSize(2) * pi/2).^2;
maskX(end-padSize(2)+1:end) = cos((0:padSize(2)-1)/padSize(2) * pi/2).^2;

W = maskY * maskX;

%% 频域坐标网格（基于填充后的尺寸）
fx = (-nCol/2 : nCol/2-1) / (nCol*PixelSize);   % cycles/m
fy = (-mRow/2 : mRow/2-1) / (mRow*PixelSize);
[u, v] = meshgrid(fx, fy);

%% 对所有帧进行倾斜位移补偿并乘以蒙版
out_img_set = cell(1,num_img);
out_img_set{1} = img_set{1} .* W;
imwrite(mat2gray(out_img_set{1}),fullfile(title_folder,"001.png"))
DistanceInterval = 0;
for i = 2:num_img
    Interval = DistanceIntervalSet(i);
    DistanceInterval = Interval + DistanceInterval;
    InputImg_fft = fftshift(fft2(img_set{i}));
    ShiftX = DistanceInterval*kx*PixelSize;
    ShiftY = DistanceInterval*ky*PixelSize;
    Phase = exp(-1i * 2*pi * ( u * ShiftX + v * ShiftY));
    out_img_set{i} = abs(ifft2(ifftshift(InputImg_fft .* Phase))) .* W;
    imwrite(mat2gray(out_img_set{i}),fullfile(title_folder,sprintf("%0.3d.png",i)))
end

MNZ_result.orig_M = M;         % 保存原始 LED X 偏移量（像素），供下游计算 FOV 覆盖蒙版
MNZ_result.orig_N = N;         % 保存原始 LED Y 偏移量（像素）
MNZ_result.orig_size = [mRow_orig, nCol_orig];  % 原始（未填充）图像尺寸
MNZ_result.M = 0;
MNZ_result.N = 0;
MNZ_result.padSize = padSize;
