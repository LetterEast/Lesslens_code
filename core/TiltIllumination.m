function [out_img_set,MNZ_result] = TiltIllumination(img_set,MNZ_result,PixelSize,DistanceIntervalSet)

M = MNZ_result.M;
N = MNZ_result.N;
Z = MNZ_result.Z;
kx = M/Z; %X方向正切值
ky = N/Z; %Y方向正切值
title_folder = "tilt_folder";
if ~exist(title_folder,"dir"),mkdir(title_folder); end

num_img = length(img_set);
[mRow, nCol] = size(img_set{1});
% --- 动态计算所需的最大位移量，以此决定边缘复制 Padding 的大小 ---
DistanceInterval = 0;
for i = 2:num_img
    DistanceInterval = DistanceInterval + DistanceIntervalSet(i);
end
max_shift_x_pixels = abs(DistanceInterval * kx);
max_shift_y_pixels = abs(DistanceInterval * ky);

% 加 20 像素余量，确保卷边(Wrap-around)只发生在 Padding 区域内
padY = ceil(max_shift_y_pixels) + 20;
padX = ceil(max_shift_x_pixels) + 20;
padY = 0;
padX = 0;

mRow_pad = mRow + 2*padY;
nCol_pad = nCol + 2*padX;

fx_pad = (-nCol_pad/2 : nCol_pad/2-1) / (nCol_pad*PixelSize);
fy_pad = (-mRow_pad/2 : mRow_pad/2-1) / (mRow_pad*PixelSize);
[u_pad, v_pad] = meshgrid(fx_pad, fy_pad);

out_img_set = cell(1,num_img);
out_img_set{1} = img_set{1};

DistanceInterval = 0;
for i = 2:num_img
    Interval = DistanceIntervalSet(i);
    DistanceInterval = Interval + DistanceInterval;

    % 1. 用 0 填充 Padding (Zero Padding)
    img_padded = padarray(img_set{i}, [padY, padX], 0, 'both');

    % 2. 在频域执行精确的亚像素平移
    InputImg_fft = fftshift(fft2(img_padded));
    ShiftX = DistanceInterval * kx * PixelSize;
    ShiftY = DistanceInterval * ky * PixelSize;
    Phase = exp(-1i * 2*pi * (u_pad * ShiftX + v_pad * ShiftY));
    shifted_padded = abs(ifft2(ifftshift(InputImg_fft .* Phase)));

    % 3. 裁剪回原图大小（卷边和衍射伪影已被安全隔离在裁剪掉的 Padding 区域）
    out_img_set{i} = shifted_padded(padY+1 : end-padY, padX+1 : end-padX);
    imwrite(mat2gray(out_img_set{i}), fullfile(title_folder, sprintf("%0.3d.png", i)));

    % 4. 计算该平面的空间有效掩膜 (Spatial Soft Mask)
    % 图像物理上向右下平移了 ShiftX / ShiftY，所以有效视场也会跟着平移。
    s_x = ShiftX / PixelSize;
    s_y = ShiftY / PixelSize;

    % ================== 额外增加四周各 5 像素的强制掩膜 (方便测试右侧伪影，测试完可注释掉) ==================
    extra_pad = 0;
    col_start = max(1 + extra_pad, round(1 + s_x) + extra_pad);
    col_end   = min(nCol - extra_pad, round(nCol + s_x) - extra_pad);
    row_start = max(1 + extra_pad, round(1 + s_y) + extra_pad);
    row_end   = min(mRow - extra_pad, round(mRow + s_y) - extra_pad);
    % =========================================================================================================

    mask_x = zeros(1, nCol);
    mask_y = zeros(mRow, 1);

    trans = 0;
    if col_end > col_start
        mask_x(col_start : col_end) = 1;
        t = min(trans, floor((col_end - col_start)/2));
        if t > 0
            ramp = 0.5 * (1 - cos(pi * (0:t-1) / t));
            if col_start > 1
                mask_x(col_start : col_start+t-1) = ramp;
            end
            if col_end < nCol
                mask_x(col_end-t+1 : col_end) = fliplr(ramp);
            end
        end
    end
    if row_end > row_start
        mask_y(row_start : row_end) = 1;
        t = min(trans, floor((row_end - row_start)/2));
        if t > 0
            ramp = 0.5 * (1 - cos(pi * (0:t-1) / t));
            if row_start > 1
                mask_y(row_start : row_start+t-1) = ramp';
            end
            if row_end < mRow
                mask_y(row_end-t+1 : row_end) = flipud(ramp');
            end
        end
    end
    MNZ_result.ValidMask{i} = mask_y * mask_x;
end


% ================== 平面1也加上额外的 5 像素强制掩膜 (测试完可注释掉恢复为全 1) ==================
MNZ_result.ValidMask{1} = zeros(mRow, nCol);
MNZ_result.ValidMask{1}(1+extra_pad : mRow-extra_pad, 1+extra_pad : nCol-extra_pad) = 1;
% MNZ_result.ValidMask{1} = ones(mRow, nCol); % <--- 原本的代码
% ==============================================================================================

MNZ_result.M = 0;
MNZ_result.N = 0;
