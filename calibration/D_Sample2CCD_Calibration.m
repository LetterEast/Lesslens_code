function [D_Sample2CCD] = D_Sample2CCD_Calibration(ImgRec,ParaD_Sample2CCD,MainPara,foldername)

%% parameter
Pre = ParaD_Sample2CCD.D_Sample2CCDPre;  %大致距离
HalfRange = ParaD_Sample2CCD.D_Sample2CCDHalfRange; %单边范围
rough = ParaD_Sample2CCD.rough; % 粗测步长

WaveLength          = MainPara.WaveLength;
PixelSize           = MainPara.PixelSize;


%% Distance Interval Calculate
D_Sample2CCD_arr = Pre + [-HalfRange:rough:HalfRange];
n_distance = length(D_Sample2CCD_arr);

%% get D_Sample2CCD
tic
ObjectPower = zeros(1,n_distance);
for i = 1 : n_distance
    Object = propGPU(ImgRec,PixelSize,WaveLength,-D_Sample2CCD_arr(:,i)); %重建的图传回Object平面
    Object_amp = abs(Object);
    Object_filt = imgaussfilt(Object_amp, 2); % 对振幅进行高斯滤波
    Object_filt = Object_amp - Object_filt + 1e-6;;
    ObjectFT = fft2(Object_filt);
    ObjectPower(:,i) = sum(sum(abs(ObjectFT)));
end

% figure, plot(D_Sample2CCD_arr,ObjectPower);
[~,index] = max(ObjectPower); % 求最大值
D_Sample2CCD = D_Sample2CCD_arr(index);

% 保存
fig = figure;
plot(D_Sample2CCD_arr,ObjectPower);
frame = getframe(fig);
img = frame2im(frame);
filefolder = fullfile(foldername,'D_sample2CCD.png');
imwrite(img, filefolder);

toc
disp(['最佳传播距离为：', num2str(D_Sample2CCD, 10), ' m']);
disp('finish D_Sample2CCD Calibration');
