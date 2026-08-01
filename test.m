clc;close all;clear all;

input_folder = 'D:\Desktop\Lesslens code\RawImgAverage';
load('D:\Desktop\data\2026\USAF_1951\2026-06-16_15-58-49\MNZ_result.mat');
output_folder = './test';
if ~exist(output_folder,'dir');mkdir(output_folder);end
M = MNZ_result.M;
N = MNZ_result.N;
Z = MNZ_result.Z;
kx = M/Z; %X方向正切值
ky = N/Z; %Y方向正切值

MNZ_result.M = 0;
MNZ_result.N = 0;

a = 0.7;b = -0.13*a + 0.6;

pad = 500;
input_img = dir(fullfile(input_folder, '*.png'));
num_img = numel(input_img);

wavelenth = 514e-9;
pixel_size = 3e-6;
Distance_step = 1e-3;

for i_img = 1:num_img
    img = im2double(imread(fullfile(input_folder, input_img(i_img).name)));
    img_pad = padarray(img,[pad,pad],0,'both');
    img_set{i_img} = img_pad;
end

% ---------------------------------------------------------
% 基于大视场的合成孔径并行迭代重建 (APRW 动量加速版)
% ---------------------------------------------------------
num_iters = 10;
Z_focus_offset = 0; % 设为 0，即默认把大画布建立在第一张传感器的 Z=0 面上，方便后续自动对焦

[H_pad, W_pad] = size(img_set{1});
H_orig = H_pad - 2*pad;
W_orig = W_pad - 2*pad;

% 使用 Tukey 窗作为探测器面的振幅替换掩膜，保护 Padding 区域保留衍射场
meas_mask = padarray(tukeywin(H_orig, 0.1) * tukeywin(W_orig, 0.1)', [pad, pad], 0, 'both');

% 计算平移量 (以像素为单位)
tx_all = round(Distance_step * (0:num_img-1) * kx);
ty_all = round(Distance_step * (0:num_img-1) * ky);

% 为了在大画布中提取子区域，平移量取反并平移到正数范围
c_offset = tx_all;
c_offset = c_offset - min(c_offset);
r_offset = ty_all;
r_offset = r_offset - min(r_offset);

% 计算合成大画布的物理总尺寸
H_large = H_pad + max(r_offset);
W_large = W_pad + max(c_offset);

% =========================================================
% 第一步：初始化大画布 (所有图像反向传播到第一面并平均)
% =========================================================
fprintf('正在初始化大视场画布...\n');
bg_val = mean(sqrt(abs(img_set{1}(:))));
Obj_Large = zeros(H_large, W_large);
Mask_Large = zeros(H_large, W_large);

for i_img = 1:num_img
    % 构造当前面对焦物体面的照明球面波
    MNZ_i = MNZ_result; MNZ_i.M = 0; MNZ_i.N = 0;
    % 因为画布在 Z_focus_offset 处，LED 到物体的距离相应减少
    MNZ_i.Z = MNZ_result.Z + (i_img-1)*Distance_step - Z_focus_offset;
    Illum_i = getMultiAngleIllum(img_set{1}, MNZ_i, pixel_size, wavelenth);

    % 原图反向传播到对焦物体面
    z_prop = (i_img-1)*Distance_step + Z_focus_offset;

    % 在传感器面上，加上一个虚拟的球面波相位以匹配光路（原代码逻辑）
    % 这里用传感器面的照明相位
    MNZ_sensor = MNZ_result; MNZ_sensor.M = 0; MNZ_sensor.N = 0;
    MNZ_sensor.Z = MNZ_result.Z + (i_img-1)*Distance_step;
    Illum_sensor = getMultiAngleIllum(img_set{1}, MNZ_sensor, pixel_size, wavelenth);

    U_backward = img_set{i_img} .* exp(1j*angle(Illum_sensor));
    if z_prop == 0
        img_back = U_backward;
    else
        img_back = propGPU(U_backward, pixel_size, wavelenth, -z_prop);
    end

    % 剥离物体面的照明相位，还原为纯物体复振幅
    SubObj_init = img_back .* exp(-1j*angle(Illum_i));

    % 拼接到大画布
    r1 = r_offset(i_img) + 1; r2 = r1 + H_pad - 1;
    c1 = c_offset(i_img) + 1; c2 = c1 + W_pad - 1;

    Obj_Large(r1:r2, c1:c2) = Obj_Large(r1:r2, c1:c2) + SubObj_init .* meas_mask;
    Mask_Large(r1:r2, c1:c2) = Mask_Large(r1:r2, c1:c2) + meas_mask;
end

% 归一化平均
Mask_Large_sum = Mask_Large;
is_zero_mask = (Mask_Large_sum == 0);
Mask_Large_sum(is_zero_mask) = 1;
Obj_Large = Obj_Large ./ Mask_Large_sum;

% 【极其关键】对于完全没有被图像覆盖的 padding 外围区域，必须赋予背景亮场 bg_val！
Obj_Large(is_zero_mask) = bg_val;

% 初始化 APRW 动量历史变量
Obj_Large_prev = Obj_Large;
Obj_Large_prev2 = Obj_Large;

% =========================================================
% 第二步：并行交替投影迭代 (APRW)
% =========================================================
fprintf('开始 APRW 并行迭代，总次数：%d\n', num_iters);

for iter = 1:num_iters
    Obj_Large_Update = zeros(H_large, W_large);

    % 并行计算每个视场的更新
    for i_img = 1:num_img
        r1 = r_offset(i_img) + 1; r2 = r1 + H_pad - 1;
        c1 = c_offset(i_img) + 1; c2 = c1 + W_pad - 1;

        SubObj = Obj_Large(r1:r2, c1:c2);

        % 构造当前面对焦物体面的照明球面波
        MNZ_i = MNZ_result; MNZ_i.M = 0; MNZ_i.N = 0;
        MNZ_i.Z = MNZ_result.Z + (i_img-1)*Distance_step - Z_focus_offset;
        Illum_i = getMultiAngleIllum(img_set{1}, MNZ_i, pixel_size, wavelenth);

        % 正向传播到传感器面
        U_forward = SubObj .* Illum_i;
        z_prop = (i_img-1)*Distance_step + Z_focus_offset;
        if z_prop == 0
            lightOnDetector = U_forward;
        else
            lightOnDetector = propGPU(U_forward, pixel_size, wavelenth, z_prop);
        end

        % 振幅替换 (加权融合)
        predictedAmp = abs(lightOnDetector);
        measuredAmp = sqrt(abs(img_set{i_img}));
        updatedAmp = meas_mask .* measuredAmp + (1 - meas_mask) .* predictedAmp;
        lightOnDetector = updatedAmp .* exp(1j * angle(lightOnDetector));

        % 反向传播
        if z_prop == 0
            U_backward = lightOnDetector;
        else
            U_backward = propGPU(lightOnDetector, pixel_size, wavelenth, -z_prop);
        end

        % 剥离照明相位，得到该视场对纯物体的预测更新
        SubObj_updated = U_backward .* exp(-1j * angle(Illum_i));

        % 累加更新量
        Obj_Large_Update(r1:r2, c1:c2) = Obj_Large_Update(r1:r2, c1:c2) + SubObj_updated .* meas_mask;
    end

    % 计算当前迭代周期的全局平均对象
    Obj_Large_Average = Obj_Large_Update ./ Mask_Large_sum;
    % 每轮迭代后同样要强制将外围填充为背景亮场，避免边缘衍射灾难
    Obj_Large_Average(is_zero_mask) = bg_val;

    % 施加 APRW 动量加速
    if iter <= 2
        % 原版 APRW 核心逻辑：前两次迭代极不稳定，禁止加动量，否则振幅会被减到负数引发灾难翻转
        Obj_Large_Next = Obj_Large_Average;
    else
        Obj_Large_Next = (1+a+b) * Obj_Large_Average - a * Obj_Large_prev - b * Obj_Large_prev2;
    end

    % 记录状态与更新
    Obj_Large_prev2 = Obj_Large_prev;
    Obj_Large_prev = Obj_Large_Next;
    Obj_Large = Obj_Large_Next;

    % 物理吸收约束 (吸收不能产生过多能量)
    Obj_Large = min(abs(Obj_Large), 1.05) .* exp(1j * angle(Obj_Large));

    fprintf('Iter %d / %d done.\n', iter, num_iters);
end

% --- 重建完成，裁剪并输出 ---
ImgRec_mean = Obj_Large(pad+1:end-pad, pad+1:end-pad);

imwrite(mat2gray(abs(ImgRec_mean)), fullfile(output_folder,'mean_reconstructed_image.png'));
figure, imshow(mat2gray(abs(ImgRec_mean))); title('Synthetic Aperture (Parallel APRW Iterative) - Z=0');
disp('done');

Object = propGPU(ImgRec_mean, pixel_size, wavelenth, -1.4e-3);
imwrite(mat2gray(abs(Object)), fullfile(output_folder,'Object_image.png'));
figure, imshow(mat2gray(abs(Object))); title('Object_image');