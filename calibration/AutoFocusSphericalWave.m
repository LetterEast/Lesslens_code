function [DistanceIntervalSet] = AutoFocusSphericalWave(img_set,numImages,MainPara)
% AutoFocusSphericalWave : 基于球面波补偿的粗‑细两阶段自动聚焦
%
% 改进：
%   1. 物理模型：动态更新 Z 并施加对应的球面波相位，补偿 CCD 移动时的曲率变化。
%   2. 评价函数：采用 FFT 频谱能量 (L1 范数) 替代 Tenengrad，极大地降低共轭像干扰。
%   3. Tracking：第 2~N 帧的搜索范围紧紧绑定在上一帧结果上，彻底杜绝位置倒挂。

%% ------------------- 参数 -------------------
output_folder = "Distance_Interval";
if ~exist(output_folder,"dir"), mkdir(output_folder); end

Z_original = MainPara.MNZ_result.Z;             % 初始 LED2CCD 距离
Pre = MainPara.ParaD_Sample2CCD.D_Sample2CCDPre; % 第1帧估计绝对距离
Interval = MainPara.ParaD_Sample2CCD.Interval;  % 预期相邻帧移动间距
ParaD_Sample2CCD = MainPara.ParaD_Sample2CCD;
PixelSize = MainPara.PixelSize;
WaveLength = MainPara.WaveLength;

rough = ParaD_Sample2CCD.rough;
if isfield(ParaD_Sample2CCD,'fine')
    fine = ParaD_Sample2CCD.fine;
else
    fine = rough/20;
end
HalfRange = ParaD_Sample2CCD.D_Sample2CCDHalfRange;

% 用于记录每一帧找到的绝对位置 (Sample2CCD)
absPos = zeros(1, numImages);

%% ------------------- 裁剪框（仅选一次） -------------------
tic
ImgRec_amp_first = sqrt(img_set{1});
figure, [~,rect] = imcrop(ImgRec_amp_first); close;

%% ------------------- 主循环 -------------------
for i_img = 1:numImages
    fprintf('\n=== Processing Image %d / %d ===\n', i_img, numImages);
    ImgRec_raw = sqrt(img_set{i_img});
    
    % --- 确定当前帧的搜索范围 ---
    if i_img == 1
        % 第 1 帧：全范围搜索
        searchCenter = Pre;
        curHalfRange = HalfRange;
    else
        % 第 2~N 帧 (Tracking)：以上一帧的真实峰值 + 预期间距 作为中心
        searchCenter = absPos(i_img-1) + Interval;
        % 极大地缩小搜索范围，容差仅设为 ± 0.05 mm（50微米）
        curHalfRange = 0.05e-3; 
    end
    
    %% ---------------- 粗搜索 ----------------
    D_coarse = searchCenter + (-curHalfRange : rough : curHalfRange);
    n_coarse = length(D_coarse);
    coarseMetric = zeros(1,n_coarse);
    for k = 1:n_coarse
        coarseMetric(k) = evalFocusMetric(D_coarse(k), ImgRec_raw, rect, Z_original, Pre, MainPara, PixelSize, WaveLength);
    end
    [~, idxCoarse] = max(coarseMetric);
    d_coarse_peak = D_coarse(idxCoarse);
    fprintf('  Coarse peak : %.6f m (%.4f mm)\n', d_coarse_peak, d_coarse_peak*1000);

    %% ---------------- 细搜索（在粗峰值 ±rough） ----------------
    D_fine = d_coarse_peak + (-rough : fine : rough);
    n_fine = length(D_fine);
    fineMetric = zeros(1,n_fine);
    for k = 1:n_fine
        fineMetric(k) = evalFocusMetric(D_fine(k), ImgRec_raw, rect, Z_original, Pre, MainPara, PixelSize, WaveLength);
    end
    [~, idxFine] = max(fineMetric);
    d_fine_peak = D_fine(idxFine);

    %% ----- 抛物线亚像素插值 -----
    if idxFine>1 && idxFine<n_fine
        y1 = fineMetric(idxFine-1);
        y2 = fineMetric(idxFine);
        y3 = fineMetric(idxFine+1);
        denom = y3 - 2*y2 + y1;
        if abs(denom) > eps
            d_best = d_fine_peak - 0.5*(y3 - y1)/denom * fine;
        else
            d_best = d_fine_peak;
        end
    else
        d_best = d_fine_peak;
    end

    absPos(i_img) = d_best;
    fprintf('  Final (interp): %.6f m (%.4f mm)\n', d_best, d_best*1000);

    %% ---------------- 保存曲线图 ----------------
    fig = figure('Visible','off');
    subplot(2,1,1);
    plot(D_coarse*1000, coarseMetric, '-b','LineWidth',1.2);
    hold on; plot(d_coarse_peak*1000, coarseMetric(idxCoarse), 'ro','MarkerSize',8,'MarkerFaceColor','r');
    title(['Image ', num2str(i_img), ' – Coarse Search (FFT Power)']);
    xlabel('Distance (mm)'); ylabel('Focus Metric'); grid on;
    subplot(2,1,2);
    plot(D_fine*1000, fineMetric, '-b','LineWidth',1.2);
    hold on; plot(d_fine_peak*1000, fineMetric(idxFine), 'r^','MarkerSize',10,'MarkerFaceColor','r');
    title(['Image ', num2str(i_img), ' – Fine Search']);
    xlabel('Distance (mm)'); ylabel('Focus Metric'); grid on;
    saveas(fig, fullfile(output_folder, ['D_sample2CCD_', num2str(i_img), '.png']));
    close(fig);
end

%% ------------------- 后处理 -------------------
% 因为使用了 Tracking，absPos 必然是严格递增的，无需强制 cummax
fprintf('\n=== Absolute Positions (mm) ===\n');
fprintf('  %.4f', absPos*1000);
fprintf('\n');

% 转换为相邻帧的间距
DistanceIntervalSet = [0, diff(absPos)];

% ------------------- 系统性校正 -------------------
if exist('calibration_factor.mat','file')
    load('calibration_factor.mat','c');
    DistanceIntervalSet = DistanceIntervalSet * c;
    fprintf('Applied calibration factor: %.6f\n', c);
end

toc
fprintf('\n=== DistanceIntervalSet (mm) ===\n');
fprintf('  %.4f', DistanceIntervalSet*1000);
fprintf('\n');

disp('finish Distance Interval Calibration');
end

%% ================= 辅助函数 ================= %%
function metric = evalFocusMetric(d_test, ImgRec_raw, rect, Z_original, Pre, MainPara, PixelSize, WaveLength)
    % 1. 动态更新 Z (因为 CCD 在移动，D_LED2CCD 会变化)
    z_test = Z_original + (d_test - Pre);
    MainPara.MNZ_result.Z = z_test;
    
    % 2. 计算当前距离对应的全局球面波
    Illum = getMultiAngleIllum(ImgRec_raw, MainPara.MNZ_result, PixelSize, WaveLength);
    
    % 3. 对全图裁剪，并乘上裁剪后的球面波相位 (抵消球面曲率)
    ImgRec_crop = imcrop(ImgRec_raw, rect);
    IllumCrop = imcrop(angle(Illum), rect);
    if ~isempty(IllumCrop)
        ImgRec = ImgRec_crop .* exp(1j .* IllumCrop);
    else
        ImgRec = ImgRec_crop;
    end
    
    % 4. 此时光场已被展平，可直接用 plane wave 模型 (propGPU) 进行反向传播
    Object = propGPU(ImgRec, PixelSize, WaveLength, -d_test);
    Amp = abs(Object);
    
    % 5. FFT 频谱能量 (减去均值滤除直流，突出高频纹理)
    % 这种度量对全息的共轭像干扰具有极强的鲁棒性
    Amp_AC = Amp - mean(Amp(:));
    ObjectFT = fft2(Amp_AC);
    metric = sum(abs(ObjectFT), 'all');
end
