function img_set = loadImg(img_folder, TargetSize)
% img_folder      : 图像所在文件夹
% n_img           : 读取的图像个数
% filename_format : 文件名格式（用于 sprintf）
% step            : 步长（可选，默认 1）
% TargetSize      : 填充后的尺寸（可选）
%                   - 标量 S  -> [S, S]
%                   - 二元向量 [Mpad, Npad]
ImageFormat = {'.png', '.jpg', '.jpeg', '.bmp', '.tif', '.tiff'};
input_folder   = img_folder;
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
n_img = length(img_folder_info);
img_set = cell(n_img, 1);



% 是否需要零填充
doPad = (nargin >= 2) && ~isempty(TargetSize);

if doPad
    if isscalar(TargetSize)
        Mpad = round(TargetSize);
        Npad = round(TargetSize);
    elseif numel(TargetSize) == 2
        Mpad = round(TargetSize(1));
        Npad = round(TargetSize(2));
    else
        error('TargetSize 必须是标量或 [Mpad Npad] 二元向量');
    end
end

for i_img = 1:n_img
    % file_index = 1+(i_img-1)*step;
    % filename   = sprintf(filename_format, file_index);

    filepath   = fullfile(img_folder_info(1).folder, img_folder_info(i_img).name);

    img = im2double(imread(filepath));

    if size(img, 3) == 3
        img = rgb2gray(img);
        fprintf("Warning! The input images are RGB.\n");
    end

    if doPad
        [m, n] = size(img);

        if Mpad < m || Npad < n
            error('TargetSize [%d %d] 不能小于原图尺寸 [%d %d].', Mpad, Npad, m, n);
        end


        canvas = zeros(Mpad, Npad, 'like', img);

        row_start = floor((Mpad - m)/2) + 1;
        col_start = floor((Npad - n)/2) + 1;
        row_end   = row_start + m - 1;
        col_end   = col_start + n - 1;
        % row_start = 1;
        % col_start =1;
        % row_end   =  m ;
        % col_end   =  n;

        canvas(row_start:row_end, col_start:col_end) = img;

        img_set{i_img} = canvas;
    else

        img_set{i_img} = img;
    end
end
end
