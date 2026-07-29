function [out_img_set,MNZ_result] = TiltIllumination(img_set,MNZ_result,PixelSize,DistanceIntervalSet)
M = MNZ_result.M;
N = MNZ_result.N;
Z = MNZ_result.Z;
kx = M/Z; %X方向正切值
ky = N/Z; %Y方向正切值
num_img = length(img_set);
[mRow, nCol] = size(img_set{1});
% --- 动态计算所需的最大位移量，以此决定大画布的大小 ---
DistanceInterval = 0;
for i = 2:num_img
    DistanceInterval = DistanceInterval + DistanceIntervalSet(i);
end
max_shift_x_pixels = abs(DistanceInterval * kx);
max_shift_y_pixels = abs(DistanceInterval * ky);
% 加 20 像素余量
padY = ceil(max_shift_y_pixels) + 200;
padX = ceil(max_shift_x_pixels) + 200;
mRow_pad = mRow + 2*padY;
nCol_pad = nCol + 2*padX;
fx_pad = (-nCol_pad/2 : nCol_pad/2-1) / (nCol_pad*PixelSize);
fy_pad = (-mRow_pad/2 : mRow_pad/2-1) / (mRow_pad*PixelSize);
[u_pad, v_pad] = meshgrid(fx_pad, fy_pad);
out_img_set = cell(1,num_img);
% 对第一张图 Padding
out_img_set{1} = padarray(img_set{1}, [padY, padX], "replicate", 'both');
% 生成第一张图的掩膜
mask_x = zeros(1, nCol_pad);
mask_y = zeros(mRow_pad, 1);
trans = 100;
col_start = padX + 1; col_end = padX + nCol;
mask_x(col_start : col_end) = 1;
t = min(trans, floor((col_end - col_start)/2));
if t > 0
    ramp = 0.5 * (1 - cos(pi * (0:t-1) / t));
    mask_x(col_start : col_start+t-1) = ramp;
    mask_x(col_end-t+1 : col_end) = fliplr(ramp);
end
row_start = padY + 1; row_end = padY + mRow;
mask_y(row_start : row_end) = 1;
t = min(trans, floor((row_end - row_start)/2));
if t > 0
    ramp = 0.5 * (1 - cos(pi * (0:t-1) / t));
    mask_y(row_start : row_start+t-1) = ramp';
    mask_y(row_end-t+1 : row_end) = flipud(ramp');
end
MNZ_result.ValidMask{1} = mask_y * mask_x;
DistanceInterval = 0;
for i = 2:num_img
    Interval = DistanceIntervalSet(i);
    DistanceInterval = Interval + DistanceInterval;

    img_padded = padarray(img_set{i}, [padY, padX], "replicate", 'both');
    % 2. 改为在空间域执行亚像素平移
    ShiftX = DistanceInterval * kx * PixelSize;
    ShiftY = DistanceInterval * ky * PixelSize;

    tx = ShiftX / PixelSize;
    ty = ShiftY / PixelSize;

    shifted_padded = imtranslate(img_padded, [tx, ty], 'bicubic', 'OutputView', 'same');
    shifted_padded = max(shifted_padded, 0); % 防止三次插值产生微小负数

    % 3. 保存
    out_img_set{i} = shifted_padded;
    % 4. 计算掩膜
    s_x = round(ShiftX / PixelSize);
    s_y = round(ShiftY / PixelSize);
    col_start = padX + 1 + s_x;
    col_end   = padX + nCol + s_x;
    row_start = padY + 1 + s_y;
    row_end   = padY + mRow + s_y;
    col_start = max(1, min(nCol_pad, col_start));
    col_end   = max(1, min(nCol_pad, col_end));
    row_start = max(1, min(mRow_pad, row_start));
    row_end   = max(1, min(mRow_pad, row_end));
    mask_x = zeros(1, nCol_pad);
    mask_y = zeros(mRow_pad, 1);
    if col_end > col_start
        mask_x(col_start : col_end) = 1;
        t = min(trans, floor((col_end - col_start)/2));
        if t > 0
            ramp = 0.5 * (1 - cos(pi * (0:t-1) / t));
            mask_x(col_start : col_start+t-1) = ramp;
            mask_x(col_end-t+1 : col_end) = fliplr(ramp);
        end
    end
    if row_end > row_start
        mask_y(row_start : row_end) = 1;
        t = min(trans, floor((row_end - row_start)/2));
        if t > 0
            ramp = 0.5 * (1 - cos(pi * (0:t-1) / t));
            mask_y(row_start : row_start+t-1) = ramp';
            mask_y(row_end-t+1 : row_end) = flipud(ramp');
        end
    end
    MNZ_result.ValidMask{i} = mask_y * mask_x;
end
MNZ_result.M = 0;
MNZ_result.N = 0;
end
