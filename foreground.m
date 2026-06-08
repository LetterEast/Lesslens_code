clc;
close all;
clear;

folder_pixel_on  = "\\192.168.2.166\d\lesslens\2026.6.8\USAF_1951_4\pixel_on\exposure_-2";
folder_pixel_off = "\\192.168.2.166\d\lesslens\2026.6.8\USAF_1951_4\pixel_off\exposure_-2";
folder_output = "\\192.168.2.166\d\lesslens\2026.6.8\USAF_1951_4\foreground_img";
if ~exist(folder_output, 'dir')
    mkdir(folder_output);
end

file_list_on = dir(fullfile(folder_pixel_on, '*.png'));
file_list_off = dir(fullfile(folder_pixel_off, '*.png'));

num_files = numel(file_list_on);
if numel(file_list_off) ~= num_files
    error('The number of files in the two folders is not equal.');
end
disp('Processing...')
for i = 1:num_files
    image_on = im2double(imread(fullfile(folder_pixel_on, file_list_on(i).name)));
    image_off = im2double(imread(fullfile(folder_pixel_off, file_list_off(i).name)));
    foreground_img =( abs(image_on - image_off));
    imwrite(foreground_img, fullfile(folder_output, file_list_on(i).name));
end
disp('Done')