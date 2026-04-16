clc; clear; close all;
addpath(genpath(pwd));

img_folder = 'D:\Desktop\data\2026\rgb\img';
MNZ_result_folder = 'D:\Desktop\data\2026\rgb\MNZ';

output_folder = 'D:\Desktop\data\2026\rgb\result';

% 过滤掉隐藏文件夹 '.' 和 '..'
dirs = dir(fullfile(img_folder, '*'));
is_valid_dir = [dirs.isdir] & ~ismember({dirs.name}, {'.','..'});
img_info = dirs(is_valid_dir);
num_folder = length(img_info);

WaveLength                   = [460e-9, 514e-9, 647e-9]; % 通常分别对应 B, G, R
MainPara.PixelSize           = 3e-6;                            % [m/pixel] sensor effective pixel size
MainPara.nIterative          = 10;                              % iteration count
MainPara.numImages           = 5;                               % number of captured positions/images
MainPara.DistanceIntervalSet = [0, ones(1, MainPara.numImages-1)] * 1e-3; % [m] distance interval between adjacent captures

MainPara.ParaD_Sample2CCD.D_Sample2CCDPre = 2e-3;        % [m]
MainPara.ParaD_Sample2CCD.D_Sample2CCDHalfRange = 1e-3;  % [m]
MainPara.ParaD_Sample2CCD.rough = 0.01e-3;               % [m] 粗测步长

iIte_record = MainPara.nIterative;
global batch_predefined_rect;
batch_predefined_rect = [];

% 用来暂存RGB三个通道的振幅/强度重建结果
color_channels = cell(1, 3);

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
    fprintf('     - 正在计算多角度照明参数...\n');
    IllumPaSet = getMultiAngleIllum(img_set{1}, MNZ_result, MainPara.PixelSize, current_wavelength);
    MainPara.IllumSet = IllumPaSet;

    fprintf('     - 正在开始 APRW 相位恢复迭代...\n');
    [Rec_nInterative, batch_predefined_rect] = getRec_APRW_batch(img_set, MainPara, MainPara.DistanceIntervalSet, iIte_record, batch_predefined_rect);

    fprintf('     - 通道 %d 重建完成。\n\n', i_folder);

    % 获取最终迭代出来的物体并转化为强度图 [0,1]
    FinalObj = Rec_nInterative;
    IntensityImg = mat2gray(abs(FinalObj));

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

% 对齐三个通道
usfac = 20;
buf1ft = fft2(color_channels{1});
for i_img = 2:numel(color_channels)
    buf2ft = fft2(color_channels{i_img});
    [output, ~] = dftregistration(buf1ft,buf2ft,usfac);
    Shift_X = output(4);
    Shift_Y = output(3);
    color_channels{i_img} = subpixelshift(color_channels{i_img}, Shift_X, Shift_Y);
end

% 如果成功跑了3个通道，将其合成为彩色RGB输出
if num_folder == 3
    fprintf('==============================================\n');
    fprintf('正在合成最终 RGB 彩色图像...\n');
    fprintf('==============================================\n');

    % 假设按 416nm(B), 514nm(G), 632nm(R) 的波长顺序存放文件夹
    color_img = cat(3, color_channels{3}, color_channels{2}, color_channels{1});

    if ~isempty(output_folder)
        if ~exist(output_folder, 'dir')
            mkdir(output_folder);
        end
        imwrite(color_img, fullfile(output_folder, 'color_fusion_result.png'));
        fprintf('     - 已成功将彩色融合结果保存至：%s\n', output_folder);
    end
    figure; imshow(color_img); title('Color Fusion Result');
    fprintf('所有任务已完成！\n');
end