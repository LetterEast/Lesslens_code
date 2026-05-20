clear; clc; close all;
addpath(genpath(pwd));
% Define base paths
image_dir    = 'D:\Desktop\data\2026\4.17data\green_usa_0.5mm\one-bg_20260417_161556';
circle_dir = 'D:\Desktop\data\2026\4.17data\green_circle_1mm\one-bg_20260417_153339';
result_base = fullfile(pwd, 'result');

if ~exist(result_base, 'dir')
    mkdir(result_base);
end

% Get list of subfolders (e.g. 1(235,395))
dirs = dir(fullfile(image_dir, '*'));
is_valid_dir = [dirs.isdir] & ~ismember({dirs.name}, {'.','..','Output_All_Reconstructions'});
valid_dirs = dirs(is_valid_dir);

global batch_predefined_rect;
batch_predefined_rect = [];

for k = 1:length(valid_dirs)
    folder_name = valid_dirs(k).name;

    RawImgFolder = fullfile(image_dir, folder_name);

    % Find MNZ_result.mat in the corresponding circle_dir
    circle_sub_dir = fullfile(circle_dir, folder_name);
    mat_files = dir(fullfile(circle_sub_dir, '**', 'MNZ_result.mat'));
    if isempty(mat_files)
        fprintf('Cannot find MNZ_result.mat in %s. Skipping...\n', circle_sub_dir);
        continue;
    end
    mnz_mat_path = fullfile(mat_files(1).folder, mat_files(1).name);

    fprintf('\n=================================\n');
    fprintf('Processing %s\n', folder_name);
    fprintf('=================================\n');

    % Parameters
    WaveLength  = 514e-9;
    PixelSize   = 3e-6;
    nIterative  = 10;
    numImages   = 5;

    iIte_record = 10;

    ParaDisInterval.DisIntervalPre =  1e-3;
    ParaDisInterval.DisIntervalHalfRange = 0e-3;
    ParaDisInterval.rough = 0.01e-3;

    ParaD_Sample2CCD.D_Sample2CCDPre = 2e-3;
    ParaD_Sample2CCD.D_Sample2CCDHalfRange = 1e-3;
    ParaD_Sample2CCD.rough = 0.01e-3;

    % 1) Raw image average function output is in OutputPath
    OutputPath = RawImgAverage(RawImgFolder, numImages);

    % 2) Load cropped images
    img_set = loadImg(OutputPath);


    % 3) Distance interval set
    a = load(mnz_mat_path);
    MNZ_result = a.MNZ_result; % load from structure

    DistanceIntervalSet = [0, ones(1,numImages-1)] * 0.5e-3;

    % 4) Calibrate Illumination
    [img_set, MNZ_result] = TiltIllumination(img_set, MNZ_result, PixelSize, DistanceIntervalSet);
    IllumPaSet = getMultiAngleIllum(img_set{1}, MNZ_result, PixelSize, WaveLength);

    % 5) Reconstruction
    MainPara.IllumSet         = IllumPaSet;
    MainPara.WaveLength       = WaveLength;
    MainPara.PixelSize        = PixelSize;
    MainPara.nIterative       = nIterative;
    MainPara.ParaD_Sample2CCD = ParaD_Sample2CCD;
    MainPara.MNZ_result       = MNZ_result;

    % Run reconstruction with predefined_rect to avoid cropping every time
    [Rec_nInterative, batch_predefined_rect] = getRec_APRW_batch(img_set, MainPara, DistanceIntervalSet, iIte_record, batch_predefined_rect);

    % Find the most recently created APRW output folder to extract the final focused object
    res_dirs = dir(fullfile('ResultFolder', 'APRW_*'));
    [~, latest_idx] = max([res_dirs.datenum]);
    latest_res_dir = fullfile(res_dirs(latest_idx).folder, res_dirs(latest_idx).name);

    iter_folder = fullfile(latest_res_dir, sprintf('Iter_%04d', nIterative));
    obj_amp_path = fullfile(iter_folder, sprintf('Object_amp_iter%04d.png', nIterative));

    if exist(obj_amp_path, 'file')
        IntensityImg = imread(obj_amp_path);
    else
        % Fallback if the object image fails to save
        FinalObj = Rec_nInterative{end};
        IntensityImg = mat2gray(abs(FinalObj));
    end

    % Save final intensity image with name 00x.png
    save_name = sprintf('%03d.png', k);
    target_res = fullfile(result_base, save_name);
    fprintf('Saving final result directly to %s\n', target_res);
    imwrite(IntensityImg, target_res);
end
fprintf('\nBatch processing completed successfully.\n');
exit;
