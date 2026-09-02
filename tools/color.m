clc; clear; close all;
addpath(genpath(pwd));

img_folder = 'D:\Desktop\data\2026\9.1\Pawpaw\img';
MNZ_result_folder = 'D:\Desktop\data\2026\9.1\Pumpkin_sec\MNZ';

output_folder = 'D:\Desktop\data\2026\9.1\Pawpaw\color_fusion_result';

% 过滤掉隐藏文件夹 '.' 和 '..'
dirs = dir(fullfile(img_folder, '*'));
is_valid_dir = [dirs.isdir] & ~ismember({dirs.name}, {'.','..'});
img_info = dirs(is_valid_dir);
num_folder = length(img_info);

WaveLength                   = [460e-9, 514e-9, 647e-9]; % 通常分别对应 B, G, R
MainPara.PixelSize           = 3e-6;                            % [m/pixel] sensor effective pixel size
MainPara.nIterative          = 6;                              % iteration count
MainPara.numImages           = 61;                               % number of captured positions/images
MainPara.DistanceIntervalSet = [0, ones(1, MainPara.numImages-1)] * 0.1e-3; % [m] distance interval between adjacent captures

MainPara.ParaD_Sample2CCD.D_Sample2CCDPre = 1.68e-3;        % [m]
MainPara.ParaD_Sample2CCD.D_Sample2CCDHalfRange = 0e-3;  % [m]
MainPara.ParaD_Sample2CCD.rough = 0.01e-3;               % [m] 粗测步长

iIte_record = MainPara.nIterative;
global batch_predefined_rect;
batch_predefined_rect = [];

% 用来暂存RGB三个通道的振幅/强度重建结果
color_channels = cell(1, 3);
all_MNZ_results = cell(1, 3);

fprintf('\n==============================================\n');
fprintf('正在启动彩色融合重建流程...\n');
fprintf('待处理波长总数：%d\n', num_folder);
fprintf('==============================================\n\n');

for i_folder = 1:num_folder
    folder_name = img_info(i_folder).name;


    folder_prefix = folder_name(1:2);
    mat_path = fullfile(MNZ_result_folder, folder_prefix, 'MNZ_result.mat');

    if ~exist(mat_path, 'file')
        error('未找到通道 %s 对应的标定文件：%s', folder_name, mat_path);
    end

    data = load(mat_path);
    MNZ_result = data.MNZ_result;


    % 波长选择
    if i_folder <= length(WaveLength)
        current_wavelength = WaveLength(i_folder);
    else
        current_wavelength = WaveLength(end);
    end
    MainPara.WaveLength = current_wavelength;

    fprintf('---> [通道 %d/%d] 正在处理文件夹：%s\n', i_folder, num_folder, folder_name);
    fprintf('     - 波长：%g nm\n', current_wavelength * 1e9);

    current_img_folder = fullfile(img_info(i_folder).folder, folder_name);

    fprintf('     - 正在加载图像...\n');
    img_set = loadImg(current_img_folder);

    fprintf('     - 正在应用倾斜照明校正...\n');
    [img_set, MNZ_result] = TiltIllumination(img_set, MNZ_result, MainPara.PixelSize, MainPara.DistanceIntervalSet);
    MainPara.MNZ_result = MNZ_result;
    all_MNZ_results{i_folder} = MNZ_result;
    fprintf('     - 正在计算多角度照明参数...\n');
    IllumPaSet = getMultiAngleIllum(img_set{1}, MNZ_result, MainPara.PixelSize, current_wavelength);
    MainPara.IllumSet = IllumPaSet;

    fprintf('     - 正在开始 APRW 相位恢复迭代...\n');
    [Rec_nInterative, batch_predefined_rect] = getRec_APRW_batch(img_set, MainPara, MainPara.DistanceIntervalSet, iIte_record, batch_predefined_rect);

    fprintf('     - 通道 %d 重建完成。\n\n', i_folder);

    % 获取最终迭代出来的 CCD 面复振幅 (包含衍射环)
    ImgRec = Rec_nInterative{end};

    fprintf('     - 正在将光场反向传播至样品面以实现对焦...\n');
    % 使用初始的样品到 CCD 距离将光场反传，消除衍射环得到清晰的物体
    D_Sample2CCD = MainPara.ParaD_Sample2CCD.D_Sample2CCDPre;
    Object_full = back2Object(ImgRec, D_Sample2CCD, MainPara);

    % 取对焦物体的振幅并转化为强度图 [0,1]
    % (由于取 abs，照明的相位项会被自动消除去，无需显式去除球面波相位)
    IntensityImg = mat2gray(abs(Object_full));

    color_channels{i_folder} = IntensityImg;
end

% % 裁剪三通道
% for i_img = 1:numel(color_channels)
%     if i_img == 1
%         [color_channels{i_img},rect] = imcrop(color_channels{i_img});
%     else
%         color_channels{i_img} = imcrop(color_channels{i_img},rect);
%     end
% end

% 对齐三个通道：根据 MNZ_result 的物理标定参数进行平移

% 2. 根据 MNZ_result 计算物理偏差并进行亚像素平移
% 必须在去除 padding 之前进行平移，这样才能利用 padding 区域重建出的信息，
% 防止平移后边缘出现黑边。
D_Sample2CCD = MainPara.ParaD_Sample2CCD.D_Sample2CCDPre;

% 注意：TiltIllumination 会将 M 和 N 置零，因此必须读取 orig_M 和 orig_N
M1 = all_MNZ_results{1}.orig_M;
N1 = all_MNZ_results{1}.orig_N;
Z1 = all_MNZ_results{1}.Z;
kx1 = M1 / Z1;
ky1 = N1 / Z1;

for i_img = 2:numel(color_channels)
    Mi = all_MNZ_results{i_img}.orig_M;
    Ni = all_MNZ_results{i_img}.orig_N;
    Zi = all_MNZ_results{i_img}.Z;
    kxi = Mi / Zi;
    kyi = Ni / Zi;

    % 计算物理平移量 (像素)
    % kx = M/Z (单位：像素/米)，所以 D_Sample2CCD(米) * kx 就是物理平移的像素数！
    % 注意符号：照明点坐标为正时，阴影向负方向偏移 (x = -D*kx)。
    % 为了将图像 i 移动到图像 1 的位置，平移量应为 (-D*kx1) - (-D*kxi) = D * (kxi - kx1)
    Shift_X = D_Sample2CCD * (kxi - kx1);
    Shift_Y = D_Sample2CCD * (kyi - ky1);

    color_channels{i_img} = subpixelshift(color_channels{i_img}, Shift_X, Shift_Y);
end

% 3. 根据有效视场掩膜 (OutputValidMask) 裁剪掉外围多余的黑边/无用区域
% 同时保证各个通道裁剪后的维度完全一致，能够通过 cat(3) 进行合并
y_min = -inf;
y_max = inf;
x_min = -inf;
x_max = inf;

for i_img = 1:numel(all_MNZ_results)
    if isfield(all_MNZ_results{i_img}, 'OutputValidMask')
        mask = all_MNZ_results{i_img}.OutputValidMask;
        [r, c] = find(mask);
        if ~isempty(r)
            padY_i = all_MNZ_results{i_img}.padSize(1);
            padX_i = all_MNZ_results{i_img}.padSize(2);

            % 映射到统一的物理坐标系 (0点在原始图像上边缘外侧)
            y_min_i = min(r) - padY_i;
            y_max_i = max(r) - padY_i;
            x_min_i = min(c) - padX_i;
            x_max_i = max(c) - padX_i;

            % 取所有通道有效区域的交集，确保都没有黑边
            y_min = max(y_min, y_min_i);
            y_max = min(y_max, y_max_i);
            x_min = max(x_min, x_min_i);
            x_max = min(x_max, x_max_i);
        end
    end
end

% 如果没有找到有效掩膜，回退到原图大小
if isinf(y_min)
    y_min = 1;
    y_max = all_MNZ_results{1}.orig_size(1);
    x_min = 1;
    x_max = all_MNZ_results{1}.orig_size(2);
end

% 考虑到平移 subpixelshift 和掩膜计算的误差，向内缩进 5 个像素，确保最边缘绝对干净
y_min = y_min + 5;
y_max = y_max - 5;
x_min = x_min + 5;
x_max = x_max - 5;

for i_img = 1:numel(color_channels)
    padY = all_MNZ_results{i_img}.padSize(1);
    padX = all_MNZ_results{i_img}.padSize(2);

    row_start = y_min + padY;
    row_end   = y_max + padY;
    col_start = x_min + padX;
    col_end   = x_max + padX;

    color_channels{i_img} = color_channels{i_img}(row_start:row_end, col_start:col_end);
end
% 如果成功跑了3个通道，将其合成为彩色RGB输出
if num_folder == 3
    fprintf('==============================================\n');
    fprintf('正在合成最终 RGB 彩色图像...\n');
    fprintf('==============================================\n');

    % 假设按 416nm(B), 514nm(G), 632nm(R) 的波长顺序存放文件夹
    color_img = cat(3, color_channels{3}, color_channels{2}, color_channels{1});

    fprintf('     - 正在应用 YUV(YCbCr) 空间颜色平均法平滑颜色噪声...\n');
    % 将 RGB 转换为 YCbCr 颜色空间 (在MATLAB中通常用YCbCr代替理论的YUV，原理相同)
    ycbcr_img = rgb2ycbcr(color_img);

    % 分离 Y (亮度), Cb, Cr (色度) 通道
    Y  = ycbcr_img(:,:,1);
    Cb = ycbcr_img(:,:,2);
    Cr = ycbcr_img(:,:,3);

    % 定义均值滤波核，窗口大小可根据噪声程度调整 (例如 5x5, 7x7)
    % 窗口越大去色差效果越强，但也可能导致颜色溢出边缘
    filter_size = 2;
    h = fspecial('average', filter_size);

    % 仅对色度通道 (Cb, Cr) 进行平均/平滑处理，保持亮度通道 (Y) 的细节不被模糊
    Cb_filtered = imfilter(Cb, h, 'replicate');
    Cr_filtered = imfilter(Cr, h, 'replicate');

    % 合并处理后的通道并转换回 RGB 空间
    ycbcr_filtered = cat(3, Y, Cb_filtered, Cr_filtered);
    color_img_yuv_filtered = ycbcr2rgb(ycbcr_filtered);

    if ~isempty(output_folder)
        if ~exist(output_folder, 'dir')
            mkdir(output_folder);
        end
        imwrite(color_img, fullfile(output_folder, 'color_fusion_result.png'));
        imwrite(color_img_yuv_filtered, fullfile(output_folder, 'color_fusion_result_yuv_filtered.png'));
        fprintf('     - 已成功将彩色融合结果保存至：%s\n', output_folder);
    end

    figure;
    subplot(1,2,1); imshow(color_img); title('原始彩色融合结果 (Original)');
    subplot(1,2,2); imshow(color_img_yuv_filtered); title('YUV颜色平均后结果 (YUV Averaged)');
    fprintf('所有任务已完成！\n');
end