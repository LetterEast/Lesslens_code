function [DistanceIntervalSet] = DisIntervalCalibration(img_set,numImages,ParaDisInterval,PixelSize,WaveLength,IllumPaSet)
%% save folder
output_folder = "Distance_Interval";
if ~exist(output_folder,"dir"),mkdir(output_folder);end
%% parameter
Pre = ParaDisInterval.DisIntervalPre;  % 大致的距离间隔
HalfRange = ParaDisInterval.DisIntervalHalfRange;  % 调整的距离范围（单边值）
rough = ParaDisInterval.rough; % 粗测步长
DistanceInterval_arr = Pre + (-HalfRange:rough:HalfRange);
nDistance = length(DistanceInterval_arr);
%% Distance Interval Calibration
tic
AmpImgSSIM = zeros(1,nDistance);
DistanceIntervalSet = zeros(1,numImages);
ImgRec = sqrt(img_set{1});
figure,[~,rect] = imcrop(ImgRec);close;

for n = 2:numImages
    % 粗调
    ImgRec = img_set{n-1}*0.5;
    ImgRec_crop = imcrop(ImgRec,rect);
    IllumPaSet_crop = imcrop(angle(IllumPaSet),rect);
    if IllumPaSet_crop
        ImgRec = ImgRec_crop .*exp(1j.*IllumPaSet_crop);
    else
        ImgRec = ImgRec_crop;
    end
   
    for i = 1:nDistance
        %% Distance Interval Calculate

        lightOnDetector   = propGPU(ImgRec,PixelSize,WaveLength,DistanceInterval_arr(i));
        simulateAmp       = abs(gather(lightOnDetector)); % 预测的振幅
                % realityAmp        = img_set{n}.^0.5;  % 实际的振幅
        realityAmp = sqrt(max(img_set{n}, 0));  % 负值置为 0 后开平方根
        real_crop = imcrop(realityAmp,rect);
        AmpImgSSIM(:,i)   = ssim(simulateAmp,real_crop);  % 两者求相似度
    end
    [~,index] = max(AmpImgSSIM);  % 求SSIM最大值
    DistanceIntervalSet(:,n) = DistanceInterval_arr(index);  % 存入距离
   
    % 保存
    fig = figure;
    plot(DistanceInterval_arr*100,AmpImgSSIM);
    frame = getframe(fig); % 获取frame
    img = frame2im(frame); % 将frame变换成imwrite函数可以识别的格式
    filename = sprintf('%d to %d.png',n,n-1);
    filefolder = fullfile(output_folder,filename);
    imwrite(img, filefolder);

    % 细调（拟合）
    %     p = polyfit(DistanceInterval_arr, AmpImgSSIM, 2);  % 距离间隔和SSIM拟合（n个点需要n-1个多项式）.但2效果才正常
    %     AmpImgSSIM_fit = polyval(p, DistanceInterval_arr_fit); % 在更精细网格上计算原始数据和多项式拟合
    %
    %     %figure,plot(DistanceInterval_arr, AmpImgSSIM, 'o', DistanceInterval_arr_fit, AmpImgSSIM_fit,'r--'); % 拟合图像
    %     %findpeaks(AmpImgSSIM_fit,DistanceInterval_arr_fit);
    %
    %     [SSIMMax,index] = max(AmpImgSSIM_fit);  % 求SSIM最大值
    %     DistanceIntervalSet(:,n) = DistanceInterval_arr_fit(index);  % 存入距离
end
toc
disp('finish Distance Interval Calibration');
