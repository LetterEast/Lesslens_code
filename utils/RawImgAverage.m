function [OutputPath] = RawImgAverage(RawImgFolder, numImages)
% RawImgAverage - 将图像按组平均，输出 RawImgSet(:,:,j)
%
% 输入：
%   RawImgFolder - 图像文件夹路径
%   numImages    - 平均图像的组数（每组均匀取 N/numImages 张图）
%   ImageFormat  - 图像后缀名（如 'png'，不带点）
%
% 输出：
%   RawImgSet(:,:,j) - 每组的平均图像（绿色通道）

    tic
    % 获取图像文件列表
    ImageFormat = {'.png', '.jpg', '.jpeg', '.bmp', '.tif', '.tiff'};
    input_folder   = RawImgFolder;
    input_folder_info = dir(input_folder);

    for num = 1 : numel(input_folder_info)
        if input_folder_info(num).isdir
            continue
        end
        [~,~,extension] = fileparts(input_folder_info(num).name);
        if any(strcmpi(extension,ImageFormat))
            ImageFormat = extension;
            break;
        end
    end
    img_folder_info = dir(fullfile(input_folder,sprintf("*%s",ImageFormat)));
    num_img          = length(img_folder_info);

    if mod(num_img, numImages) ~= 0
        warning('图像总数不能被组数整除，将忽略多余图像');
    end

    numPerGroup = floor(num_img / numImages);

    % 获取图像尺寸
    sampleImg = im2double(imread(fullfile(RawImgFolder, img_folder_info(1).name)));
    if size(sampleImg, 3) == 1
        [H, W] = size(sampleImg);
    else
        [H, W] = size(sampleImg(:,:,2));
    end

    RawImgSet = zeros(H, W, numImages);

    % ======= 分组平均 =======
    for j = 1:numImages
        fprintf('处理第 %d/%d 组...\n', j, numImages);
        sumImg = zeros(H, W);

        for k = 1:numPerGroup
            idx = (j-1)*numPerGroup + k;
            img = im2double(imread(fullfile(RawImgFolder, img_folder_info(idx).name)));

            if size(img, 3) == 1
                green = img;
            else
                green = img(:,:,2);  % 提取绿色通道
            end

            sumImg = sumImg + green;
        end

        RawImgSet(:,:,j) = sumImg / numPerGroup;
    end

    % ======= 保存结果 =======
    OutputPath = '.\RawImgAverage';
    if ~exist(OutputPath, 'dir')
        mkdir(OutputPath);
    end

    for g = 1:numImages
        filename = fullfile(OutputPath, sprintf('%03d%s', g, ImageFormat));
        imwrite(RawImgSet(:,:,g), filename);
    end

    toc
    disp("完成图像分组平均和保存。");
end
