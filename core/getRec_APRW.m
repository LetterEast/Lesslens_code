function [Rec_nInterative] = getRec_APRW(img_set, MainPara, DistanceIntervalSet, iIte_record)
%GETREC_APRW  Multi-distance phase retrieval with weighted feedback (APRW/MDPRF style)
%
% Core idea:
%   For each iteration k, forward propagate current object estimate to each
%   measurement plane, replace amplitude with measured sqrt(I), back propagate
%   to object plane to obtain guesses g_k(n), then apply weighted feedback:
%
%   g~_k(n) = (1 + a + b) * g_k(n) - a * g~_{k-1}(n) - b * g~_{k-2}(n)
%   Feedback starts when k > 2 (k=1,2 use g~_k(n)=g_k(n)).
%   (This matches the paper’s MDPRF weighted feedback description.)
%
% Inputs:
%   img_set            : cell, measured intensity images {1..N}, size MxN (double)
%   MainPara           : struct with fields:
%                        IllumSet, WaveLength, PixelSize
%                        nIterative, ParaD_Sample2CCD
%   DistanceIntervalSet: vector, z-step list (meters). Typically:
%                        DistanceIntervalSet(1)=0, and others are step sizes.
%                        The propagation distance to plane i is sum(DistanceIntervalSet(1:i)).
%   iIte_record        : record period (iterations). If omitted, record last iter only.
%
% Output:
%   Rec_nInterative    : cell, recorded complex field ImgRec at each record point

%% ---------- Basic checks & record setting ----------
nIterative = MainPara.nIterative;

if nargin < 4 || isempty(iIte_record)
    iIte_record = nIterative; % default: only record last iteration
end

record_num = floor(nIterative / iIte_record);
Rec_nInterative = cell(record_num, 1);
count = 1;

%% ---------- Parameters setting ----------
IllumPattern     = MainPara.IllumSet;
WaveLength       = MainPara.WaveLength;   % [m]
PixelSize        = MainPara.PixelSize;    % [m/pixel]
ParaD_Sample2CCD = MainPara.ParaD_Sample2CCD;  % (kept for future calibration)
MNZ_result       = MainPara.MNZ_result;

D_LED2CCD = MNZ_result.Z;

n_img = length(img_set);
[mRow, nCol] = size(img_set{1});

% If you padded/canvas-embedded images before propagation, MarginX/Y can crop back
MarginX = 0;
MarginY = 0;
use_crop = true;

%% ---------- Weighted feedback coefficients ----------

a = 0.7;
b = -0.13*a + 0.6;   %

%% ---------- CCTV (Complex Total Variation) parameters ----------
use_tv        = true;   % Toggle TV denoising
lam_tv        = 1e-2;   % Regularization parameter (tune this, e.g., 1e-3 to 5e-2)
gam_tv        = 2;      % Gradient step size
n_subiters_tv = 20;     % TV inner loops

%% ---------- Support Constraint Parameters ----------
use_support = true;
support_pad = 30;                                  % 边缘有多少个像素受到支撑域压制
win_r       = tukeywin(mRow, 2*support_pad/mRow);
win_c       = tukeywin(nCol, 2*support_pad/nCol);
SupportMask = win_r * win_c';
%% ---------- Precompute propagation distances to each plane ----------
% zSet(i) = distance from planei
DistanceIntervalSet = DistanceIntervalSet(:).';
if length(DistanceIntervalSet) < n_img
    error('DistanceIntervalSet length (%d) must be >= number of images (%d).', length(DistanceIntervalSet), n_img);
end

zSet = zeros(n_img, 1);
DistanceInterval = 0;
for iImg = 2:n_img
    DistanceInterval = DistanceInterval + DistanceIntervalSet(iImg);
    zSet(iImg) = DistanceInterval;
end

%% ---------- Create a unique result root folder ----------
runStamp = datestr(now, 'yyyymmdd_HHMMSS');
Parent_foldername = fullfile(".\ResultFolder", sprintf("APRW_%s_N%d_I%d_rec%d", runStamp, n_img, nIterative, iIte_record));
if ~exist(Parent_foldername, 'dir'); mkdir(Parent_foldername); end

% ---------- Initialization ----------
% Initialize complex field estimate at plane 1 using sqrt of first intensity
ImgRec = sqrt(img_set{1});
ImgRec = ImgRec .* IllumPattern;  % apply illumination modulation if needed

% Store modulated guesses from previous iterations for weighted feedback
ImgRec_record1 = cell(1, n_img);  % g~_{k-1}(n)
ImgRec_record2 = cell(1, n_img);  % g~_{k-2}(n)


%% ---------- Visualization ----------
hFig = figure('Name', 'APRW / Weighted Feedback Monitor');
ax1 = subplot(1,2,1, 'Parent', hFig);
ax2 = subplot(1,2,2, 'Parent', hFig);

% placeholders
hAmp = imshow(zeros(mRow,nCol), [], 'Parent', ax1); title(ax1, 'Amplitude');
hPhs = imshow(zeros(mRow,nCol), [], 'Parent', ax2); title(ax2, 'Phase');

%% ---------- Precompute Illum_CCD (True Spherical Wave at CCD plane) ----------
k0_wave_ccd = 2 * pi / WaveLength;
x_vec_ccd = -nCol/2 * PixelSize : PixelSize : (nCol/2 - 1) * PixelSize;
y_vec_ccd = -mRow/2 * PixelSize : PixelSize : (mRow/2 - 1) * PixelSize;
[X_ccd, Y_ccd] = meshgrid(x_vec_ccd, y_vec_ccd);
LEDX_ccd = MNZ_result.M * PixelSize;
LEDY_ccd = MNZ_result.N * PixelSize;
r_ccd = sqrt((X_ccd - LEDX_ccd).^2 + (Y_ccd - LEDY_ccd).^2 + D_LED2CCD^2);
Illum_CCD = exp(1j * k0_wave_ccd .* r_ccd);

%% ---------- Reconstruction process ----------
tic;
DistanceInterval = 0;
zSet = zeros(n_img,1);
R_history = zeros(nIterative, 1);  % 用于记录每一次迭代的 R-factor

for iIte = 1 : nIterative
    % One iteration: build object guesses for all planes 1..N
    ImgRec_record = cell(1, n_img);

    err_num = 0;
    err_den = 0;

    for iImg = 1 : n_img
        if iIte == 1
            Interval         = DistanceIntervalSet(iImg);
            DistanceInterval = Interval + DistanceInterval;
            zSet(iImg)       = DistanceInterval;

        end


        % Current propagation distance from plane#1 to plane#iImg
        z    = zSet(iImg);  % [m]
        z_mm = z * 1e3;     % [mm] for display

        fprintf('Plane 1 -> %d propagation distance: %0.3f mm\n', iImg, z_mm);

        % (1) Forward propagate to detector plane iImg
        x = -MNZ_result.M*PixelSize*(z/D_LED2CCD);
        y = -MNZ_result.N*PixelSize*(z/D_LED2CCD);

        lightOnDetector = propGPU(ImgRec, PixelSize, WaveLength, z, x, y);

        % (2) Amplitude replacement
        thisAmpImg      = img_set{iImg}.^0.5;

        % --- 计算 R-factor 累加值 (使用 || A_hat - A ||_F^2 / || A ||_F^2) ---
        predictedAmp = abs(lightOnDetector);
        thisAmpImg_gpu = cast(thisAmpImg, 'like', predictedAmp);
        diffAmp = predictedAmp - thisAmpImg_gpu;
        err_num = err_num + gather(sum(diffAmp(:).^2));
        err_den = err_den + gather(sum(thisAmpImg_gpu(:).^2));

        lightOnDetector = thisAmpImg_gpu .* exp(1j .* angle(lightOnDetector));

        % (3) Back propagate to object plane => g_k(n)
        ImgRec_record{iImg} = propGPU(lightOnDetector, PixelSize, WaveLength, -z, -x, -y);

        % (4) Weighted feedback (start when k > 2)
        if iIte <= 2
            g = ImgRec_record{iImg};  % no feedback for first two iterations
        else
            g = (1 + a + b) .* ImgRec_record{iImg} ...
                - a .* ImgRec_record1{iImg} ...
                - b .* ImgRec_record2{iImg};
        end
        ImgRec_record{iImg} = g;

        % ----- Display amplitude & phase for current plane guess -----
        roi = g(MarginY+1:MarginY+mRow, MarginX+1:MarginX+nCol);

        AmpImgRec = gather(abs(roi));
        PhsImgRec = gather(angle(roi));

        AmpImgRec  = mat2gray(AmpImgRec);
        PhsImgRec  = mat2gray(PhsImgRec);
        title_name = sprintf('APRW (weighted feedback), Iter=%d, Img=%d, z=%0.3f mm', iIte, iImg, z_mm);

        set(hAmp, 'CData', AmpImgRec); xlabel(ax1, 'Amplitude'); title(ax1, title_name);
        set(hPhs, 'CData', PhsImgRec); xlabel(ax2, 'Phase');

        drawnow;
    end

    % Update history buffers for weighted feedback
    if iIte == 1
        ImgRec_record1 = ImgRec_record;
    else
        ImgRec_record2 = ImgRec_record1;
        ImgRec_record1 = ImgRec_record;
    end
    Stack  = cat(3, ImgRec_record1{:});
    % ImgRec = mean(Stack, 3);
    bg_val = mean(img_set{1}(:).^0.5);

    % 使用掩膜进行加权平均 (Mask-weighted Feedback)
    MaskStack                = cat(3, MNZ_result.ValidMask{:});
    TotalWeight              = sum(MaskStack, 3);
    ImgRec                   = (sum(Stack .* MaskStack, 3) + 1e-2 * Illum_CCD) ./ (TotalWeight + 1e-2);
    ImgRec(TotalWeight == 0) = bg_val;

    % ================== 物理约束 (Physical Constraints) ==================
    % Absorption Constraint (吸收约束):
    % 真实的纯相位或弱吸收物体，其透射率振幅物理上不能超过 1（不能放大光）。
    % 给予 1.05 的微小宽容度，防止算法为了强行拟合误差而在图像中产生极其刺眼的高亮噪点。
    ImgRec = min(abs(ImgRec), 1.05) .* exp(1j * angle(ImgRec));

    % Support Constraint (支撑域约束):
    if use_support
        ImgRec = ImgRec .* cast(SupportMask, 'like', ImgRec);
    end

    % --- 评价当前迭代的 R-factor ---
    R_history(iIte) = sqrt(err_num / err_den);
    fprintf('=== Iteration %d: R-factor = %.6f ===\n', iIte, R_history(iIte));

    % --- 判断是否收敛 ---
    is_converged = false;
    if iIte > 2
        d1 = abs(R_history(iIte) - R_history(iIte-1));
        d2 = abs(R_history(iIte-1) - R_history(iIte-2));
        % 如果连续两次迭代 R-factor 的变化极小，认为收敛
        if d1 < 5e-3 && d2 < 5e-3
            is_converged = true;
            fprintf('R-factor 变化极小，算法收敛，将执行最后一次保存并提前停止！\n');
        end
    end

    %% ---------- Record & save results ----------
    % 只要满足记录条件，或者是最后一次迭代，或者是已经收敛，就执行保存
    if (mod(iIte, iIte_record) == 0) || (iIte == nIterative) || is_converged

        % Subfolder for this recorded iteration
        foldername = fullfile(Parent_foldername, sprintf("Iter_%04d", iIte));
        if ~exist(foldername, 'dir'); mkdir(foldername); end

        % Store complex field (at reconstruction plane) in output cell
        Rec_nInterative{count} = ImgRec;

        % calibrated sample-to-CCD distance

        if iIte == iIte_record & use_crop
            Ishow = gather(mat2gray(abs(ImgRec)));
            figure,[~, rect] = imcrop(Ishow);close
        end
        if use_crop
            rect = round(rect);                 % [x y w h]
            x1 = max(1, floor(rect(1)) + 1);
            y1 = max(1, floor(rect(2)) + 1);
            w  = max(1, rect(3));
            h  = max(1, rect(4));
            x2 = min(size(ImgRec,2), x1 + w - 1);
            y2 = min(size(ImgRec,1), y1 + h - 1);
            ImgRec_crop = ImgRec(y1:y2, x1:x2);
        else
            ImgRec_crop = ImgRec;
        end

        D_Sample2CCD = D_Sample2CCD_Calibration(ImgRec_crop, ParaD_Sample2CCD, MainPara, foldername);


        % Save metadata for reproducibility
        metaPath = fullfile(foldername, "meta.mat");
        save(metaPath, 'MainPara', 'DistanceIntervalSet', 'zSet', 'a', 'b', 'D_Sample2CCD', 'iIte');

        % Back to object plane and save Object

        save(fullfile(foldername,"ImgRec.mat"), 'ImgRec');
        % 1.反向传播到样品面
        Object_full = back2Object(ImgRec, D_Sample2CCD, MainPara);

        % 2. 计算在“样品面”处的真实照明球面波
        % 样品面距离点光源 (LED) 的距离为 D_LED2CCD - D_Sample2CCD
        D_LED2Sample         = D_LED2CCD - D_Sample2CCD;
        k0_wave              = 2*pi/WaveLength;
        x_vec                = -nCol/2*(PixelSize):PixelSize:(nCol/2-1)*PixelSize;
        y_vec                = -mRow/2*PixelSize:PixelSize:(mRow/2-1)*PixelSize;
        [Object_X, Object_Y] = meshgrid(x_vec, y_vec);
        LEDX                 = MNZ_result.M * PixelSize;
        LEDY                 = MNZ_result.N * PixelSize;

        r_sample     = sqrt((Object_X - LEDX).^2 + (Object_Y - LEDY).^2 + D_LED2Sample^2);
        Illum_Sample = exp(1j * k0_wave .* r_sample);

        % 3. 在样品面上，去除球面波的相位，得到真实的纯物体透射率
        Object = Object_full .* exp(-1j * angle(Illum_Sample));


        if use_tv
            v_est = zeros(size(Object,1), size(Object,2), 2, 'like', Object);
            w_est = zeros(size(Object,1), size(Object,2), 2, 'like', Object);

            for subiter = 1:n_subiters_tv
                w_next = v_est + (1 / (8 * gam_tv)) * Df( Object  - gam_tv * DTf(v_est) );
                w_next = min(abs(w_next), lam_tv) .* exp(1j * angle(w_next));
                v_est  = w_next + (subiter / (subiter + 3)) * (w_next - w_est);
                w_est  = w_next;
            end
            Object = Object - gam_tv * DTf(w_est);
            Object = min(abs(Object), 1.05) .* exp(1j * angle(Object));
        end
        if use_support
            ImgRec = ImgRec .* cast(SupportMask, 'like', ImgRec);
        end

        save(fullfile(foldername, sprintf("Object_iter%04d.mat", iIte)), 'Object');


        % 3. 外围盲区相位置零
        Object(TotalWeight == 0) = abs(Object(TotalWeight == 0));
        % ---------------------------------------------

        % Save amplitude / phase images
        ampObj = abs(Object);
        phsObj = angle(Object);


        % 使用 imadjust 拉伸对比度
        amp_gray = (mat2gray(ampObj));
        phs_uint8 = (mat2gray(phsObj));


        pad_r      = max(1, round(size(phs_uint8, 1) * 0.25));
        pad_c      = max(1, round(size(phs_uint8, 2) * 0.25));
        phs_center = phs_uint8(pad_r:end-pad_r, pad_c:end-pad_c);

        v_min  = prctile(double(phs_center(:)), 1);
        v_max  = prctile(double(phs_center(:)), 99);
        if v_max == v_min
            v_max  = v_min + 1e-5;
        end

        pad_r      = max(1, round(size(amp_gray, 1) * 0.2));
        pad_c      = max(1, round(size(amp_gray, 2) * 0.2));
        amp_center = amp_gray(pad_r:end-pad_r, pad_c:end-pad_c);

        amp_min = prctile(double(amp_center(:)), 1);
        amp_max = prctile(double(amp_center(:)), 99);
        if amp_max == amp_min
            amp_max = amp_min + 1e-5;
        end

        % 归一化并截断
        data_norm = (double(phs_uint8) - v_min) / (v_max - v_min);
        data_norm(data_norm < 0) = 0;
        data_norm(data_norm > 1) = 1;

        amp_norm = (double(amp_gray) - amp_min) / (amp_max - amp_min);
        amp_norm(amp_norm < 0) = 0;
        amp_norm(amp_norm > 1) = 1;

        phs_gray_stretched = mat2gray(data_norm);
        amp_gray_stretched = mat2gray(amp_norm);


        imwrite(amp_gray_stretched, fullfile(foldername, sprintf("Object_amp_iter%04d.png", iIte)));
        imwrite(phs_gray_stretched, fullfile(foldername, sprintf("Object_phs_iter%04d.png", iIte)));

        % Save phase with hot colormap
        cmap = hot(256);
        img_rgb = ind2rgb(round(data_norm * 255), cmap);
        imwrite(img_rgb, fullfile(foldername, sprintf("Object_phs_hot_iter%04d.png", iIte)));

        % Quick view
        hFigObj = figure('Name', sprintf('Object @ Iter %d', iIte));
        subplot(1,3,1), imshow(amp_gray, []), xlabel("Object_amp"), axis image
        subplot(1,3,2), imshow(phs_gray_stretched), xlabel("Object_phs"), axis image
        subplot(1,3,3), imshow(img_rgb), xlabel("Object_phs (hot)"), axis image
        drawnow;
        close(hFigObj);

        count = count + 1;
    end

    % --- 如果判定已收敛，跳出迭代循环 ---
    if exist('is_converged', 'var') && is_converged
        break;
    end
end

% 绘制收敛曲线
figure('Name', 'R-factor Convergence');
plot(1:iIte, R_history(1:iIte), '-o', 'LineWidth', 2);
title('R-factor Convergence Curve');
xlabel('Iteration'); ylabel('R-factor (L2 Norm)'); grid on;

toc;
disp('finish reconstruction');
end

%% ---------- TV Denoising Helper Functions ----------
function w = Df(x)
% 2D forward difference
w = cat(3, x(1:end,:) - x([2:end,end],:), x(:,1:end) - x(:,[2:end,end]));
end

function u = DTf(w)
% Transpose of the 2D difference
u1 = w(:,:,1) - w([end,1:end-1],:,1);
u1(1,:) = w(1,:,1);
u1(end,:) = -w(end-1,:,1);

u2 = w(:,:,2) - w(:,[end,1:end-1],2);
u2(:,1) = w(:,1,2);
u2(:,end) = -w(:,end-1,2);

u = u1 + u2;
end
