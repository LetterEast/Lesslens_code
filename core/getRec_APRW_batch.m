function [Rec_nInterative, rect] = getRec_APRW_batch(img_set, MainPara, DistanceIntervalSet, iIte_record, predefined_rect)
%GETREC_APRW_BATCH  Multi-distance phase retrieval with weighted feedback (APRW/MDPRF style)
%
% Core idea:
%   For each iteration k, forward propagate current object estimate to each
%   measurement plane, replace amplitude with measured sqrt(I), back propagate
%   to object plane to obtain guesses g_k(n), then apply weighted feedback:
%
%   g~_k(n) = (1 + a + b) * g_k(n) - a * g~_{k-1}(n) - b * g~_{k-2}(n)
%   Feedback starts when k > 2 (k=1,2 use g~_k(n)=g_k(n)).
%   (This matches the paper's MDPRF weighted feedback description.)
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
%   predefined_rect    : optional predefined crop rect
%
% Output:
%   Rec_nInterative    : cell, recorded complex field ImgRec at each record point
%   rect               : the crop rect used

%% ---------- Basic checks & record setting ----------
nIterative = MainPara.nIterative;

if nargin < 4 || isempty(iIte_record)
    iIte_record = nIterative; % default: only record last iteration
end

record_num = floor(nIterative / iIte_record);
Rec_nInterative = cell(record_num, 1);
count = 1;
rect = [];
if nargin < 5
    predefined_rect = [];
end

%% ---------- Parameters setting ----------
IllumPattern     = MainPara.IllumSet;
WaveLength       = MainPara.WaveLength;   % [m]
PixelSize        = MainPara.PixelSize;    % [m/pixel]
ParaD_Sample2CCD = MainPara.ParaD_Sample2CCD;  % (kept for future calibration)
MNZ_result       = MainPara.MNZ_result;

D_LED2CCD = MNZ_result.Z;

n_img = length(img_set);
[mRow, nCol] = size(img_set{1});
propValidSize = [mRow, nCol];
if isfield(MNZ_result, 'orig_size')
    propValidSize = MNZ_result.orig_size;
end

% If you padded/canvas-embedded images before propagation, MarginX/Y can crop back
MarginX = 0;
MarginY = 0;
use_crop = true;
outputOptions = struct( ...
    'cropToValidFOV', true, ...
    'validMaskThreshold', 0.01, ...
    'zeroFillInvalid', true);
if isfield(MainPara, 'Output')
    optionNames = fieldnames(MainPara.Output);
    for iOption = 1:numel(optionNames)
        outputOptions.(optionNames{iOption}) = MainPara.Output.(optionNames{iOption});
    end
end

%% ---------- Weighted feedback coefficients ----------

a = 0.7;
b = -0.13*a + 0.6;   %
% a = 0;b=0;  % Disable weighted feedback for testing
%% ---------- CCTV (Complex Total Variation) parameters ----------
use_tv        = false;   % Toggle TV denoising
lam_tv        = 1e-3;   % Regularization parameter (tune this, e.g., 1e-3 to 5e-2)
gam_tv        = 2;      % Gradient step size
n_subiters_tv = 20;     % TV inner loops

%% ---------- Support Constraint Parameters ----------
use_support = false;
% SUPPORT-WINDOW-FIX-OLD: support_pad = 30; % Original line had win_r appended after the comment.
% >>> SUPPORT-WINDOW-FIX-BEGIN (Codex)
support_pad = 30;
win_r       = tukeywin(mRow, 2*support_pad/mRow);
% <<< SUPPORT-WINDOW-FIX-END (Codex)
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

% >>> FULL-CANVAS-INITIALIZATION-OLD-BEGIN
% % ---------- Initialization ----------
% % Initialize complex field estimate at plane 1 using sqrt of first intensity
% ImgRec = sqrt(img_set{1});
% ImgRec = ImgRec .* IllumPattern;  % apply illumination modulation if needed
% <<< FULL-CANVAS-INITIALIZATION-OLD-END

% >>> MASKED-INITIALIZATION-TEST-BEGIN (Codex)
initialAmp = sqrt(img_set{1});
% >>> SOFT-INITIALIZATION-MASK-OLD-BEGIN
% if isfield(MNZ_result, 'ValidMask') && ~isempty(MNZ_result.ValidMask)
%     initialMask = cast(MNZ_result.ValidMask{1}, 'like', initialAmp);
% else
%     initialMask = ones(size(initialAmp), 'like', initialAmp);
% end
% <<< SOFT-INITIALIZATION-MASK-OLD-END
% >>> HARD-INITIALIZATION-MASK-TEST-BEGIN (Codex)
if isfield(MNZ_result, 'ValidMaskHard') && ~isempty(MNZ_result.ValidMaskHard)
    initialMask = cast(MNZ_result.ValidMaskHard{1}, 'like', initialAmp);
else
    initialMask = ones(size(initialAmp), 'like', initialAmp);
end
% <<< HARD-INITIALIZATION-MASK-TEST-END (Codex)
initialMask = min(max(initialMask, 0), 1);
initialBackgroundAmp = sum(initialMask(:) .* initialAmp(:)) ./ ...
    max(sum(initialMask(:)), eps('like', initialAmp));
initialAmp = initialMask .* initialAmp + ...
    (1 - initialMask) .* initialBackgroundAmp;
ImgRec = initialAmp .* IllumPattern;
% <<< MASKED-INITIALIZATION-TEST-END (Codex)

% Store modulated guesses from previous iterations for weighted feedback
ImgRec_record1 = cell(1, n_img);  % g~_{k-1}(n)
ImgRec_record2 = cell(1, n_img);  % g~_{k-2}(n)


%% ---------- Visualization ----------
% hFig = figure('Name', 'APRW / Weighted Feedback Monitor');
% ax1 = subplot(1,2,1, 'Parent', hFig);
% ax2 = subplot(1,2,2, 'Parent', hFig);
% 
% % placeholders
% hAmp = imshow(zeros(mRow,nCol), [], 'Parent', ax1); title(ax1, 'Amplitude');
% hPhs = imshow(zeros(mRow,nCol), [], 'Parent', ax2); title(ax2, 'Phase');

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

        lightOnDetector = propGPU( ...
            ImgRec, PixelSize, WaveLength, z, x, y, propValidSize);

        % (2) Amplitude replacement
        thisAmpImg      = img_set{iImg}.^0.5;



        % >>> MASKED-AMPLITUDE-REPLACEMENT-TEST-BEGIN (Codex)
        predictedAmp = abs(lightOnDetector);
        thisAmpImg_gpu = cast(thisAmpImg, 'like', predictedAmp);

        if isfield(MNZ_result, 'ValidMaskHard') && ...
                numel(MNZ_result.ValidMaskHard) >= iImg
            measurementMask = cast(MNZ_result.ValidMaskHard{iImg}, ...
                'like', predictedAmp);
        else
            measurementMask = ones(size(predictedAmp), 'like', predictedAmp);
        end
        % <<< HARD-AMPLITUDE-MASK-TEST-END (Codex)
        measurementMask = min(max(measurementMask, 0), 1);

        % Evaluate R-factor only where measured information is valid.
        diffAmp = predictedAmp - thisAmpImg_gpu;
        err_num = err_num + gather(sum( ...
            measurementMask(:) .* abs(diffAmp(:)).^2));
        err_den = err_den + gather(sum( ...
            measurementMask(:) .* abs(thisAmpImg_gpu(:)).^2));

        % Valid area: measured amplitude. Replicated padding: predicted
        % amplitude. The cosine mask provides a soft transition in between.
        updatedAmp = measurementMask .* thisAmpImg_gpu + ...
            (1 - measurementMask) .* predictedAmp;
        lightOnDetector = updatedAmp .* exp(1j .* angle(lightOnDetector));
        % <<< MASKED-AMPLITUDE-REPLACEMENT-TEST-END (Codex)

        % (3) Back propagate to object plane => g_k(n)
        ImgRec_record{iImg} = propGPU( ...
            lightOnDetector, PixelSize, WaveLength, -z, -x, -y, propValidSize);

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
%         roi = g(MarginY+1:MarginY+mRow, MarginX+1:MarginX+nCol);
% 
%         AmpImgRec = gather(abs(roi));
%         PhsImgRec = gather(angle(roi));
% 
%         AmpImgRec  = mat2gray(AmpImgRec);
%         PhsImgRec  = mat2gray(PhsImgRec);
%         title_name = sprintf('APRW (weighted feedback), Iter=%d, Img=%d, z=%0.3f mm', iIte, iImg, z_mm);
% 
%         set(hAmp, 'CData', AmpImgRec); xlabel(ax1, 'Amplitude'); title(ax1, title_name);
%         set(hPhs, 'CData', PhsImgRec); xlabel(ax2, 'Phase');
% 
%         drawnow;
    end

    % Update history buffers for weighted feedback
    if iIte == 1
        ImgRec_record1 = ImgRec_record;
    else
        ImgRec_record2 = ImgRec_record1;
        ImgRec_record1 = ImgRec_record;
    end
    Stack  = cat(3, ImgRec_record1{:});
    ImgRec = mean(Stack, 3);
    % bg_val = mean(img_set{1}(:).^0.5);

    % (Mask-weighted Feedback)
    MaskStack                = cat(3, MNZ_result.ValidMask{:});
    TotalWeight              = sum(MaskStack, 3);


    % ================== 物理约束 (Physical Constraints) ==================
    % Absorption Constraint (吸收约束):
    % 真实的纯相位或弱吸收物体，其透射率振幅物理上不能超过 1（不能放大光）。
    % 给予 1.05 的微小宽容度，防止算法为了强行拟合误差而在图像中产生极其刺眼的高亮噪点。
    % >>> IMGREC-CONSTRAINT-FIX-BEGIN (Codex)
    ImgRec = min(abs(ImgRec), 1.05) .* exp(1j * angle(ImgRec));
    % <<< IMGREC-CONSTRAINT-FIX-END (Codex)

    % Support Constraint
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
        if is_converged || (iIte == nIterative)
            save_iter = nIterative; % Always save final as Iter_0010 so batch_run can find it
        else
            save_iter = iIte;
        end
        foldername = fullfile(Parent_foldername, sprintf("Iter_%04d", save_iter));
        if ~exist(foldername, 'dir'); mkdir(foldername); end

        % Store complex field (at reconstruction plane) in output cell
        Rec_nInterative{count} = ImgRec;

        % calibrated sample-to-CCD distance

        if use_crop && isempty(rect)
            if isempty(predefined_rect)
                % In batch mode, if no predefined rect is given, crop to center 80% to avoid manual popup
                rect = [round(size(ImgRec,2)*0.1), round(size(ImgRec,1)*0.1), round(size(ImgRec,2)*0.8), round(size(ImgRec,1)*0.8)];
            else
                rect = predefined_rect;
            end
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
        % >>> OBJECT-INIT-FIX-BEGIN (Codex)
        Object = Object_full .* exp(-1j * angle(Illum_Sample));
        % <<< OBJECT-INIT-FIX-END (Codex)

        % Save phase at three stages with one fixed [-pi, pi] display scale.
        % This distinguishes iterative-boundary artifacts from back-propagation
        % artifacts and illumination-removal mismatch.
        phase_ImgRec      = gather(angle(ImgRec));
        phase_ObjectFull  = gather(angle(Object_full));
        phase_Object      = gather(angle(Object));
        save(fullfile(foldername, sprintf( ...
            'PhaseStages_iter%04d.mat', iIte)), ...
            'phase_ImgRec', 'phase_ObjectFull', 'phase_Object');

        phaseStageNames = {'ImgRec_CCD', 'ObjectFull_beforeIllumRemoval', ...
            'Object_afterIllumRemoval'};
        phaseStageData = {phase_ImgRec, phase_ObjectFull, phase_Object};
        for iPhaseStage = 1:numel(phaseStageData)
            phaseFixed = uint16(max(0, min(1, ...
                (phaseStageData{iPhaseStage} + pi) / (2*pi))) * 65535);
            imwrite(phaseFixed, fullfile(foldername, sprintf( ...
                'Phase_%s_iter%04d.png', phaseStageNames{iPhaseStage}, iIte)));
        end

        if use_tv
            if isfield(MNZ_result, 'ValidMaskHard')
                coverageCount = sum(cat(3, MNZ_result.ValidMaskHard{:}), 3);
            else
                coverageCount = double(TotalWeight > outputOptions.validMaskThreshold);
            end
            maxCoverage = max(coverageCount(:));
            confidenceMap = coverageCount ./ max(maxCoverage, 1);
            reliableTVMask = coverageCount > 0;

            % Multi-plane disagreement is retained, but it is gated by the
            % crop-boundary risk and therefore cannot strengthen TV at center.
            residualFloor = 0.01 * mean(abs(ImgRec(reliableTVMask)).^2, 'all');
            residualDenominator = abs(ImgRec).^2 + residualFloor + eps('like', real(ImgRec));
            residualEnergy = sum(MaskStack .* abs(Stack - ImgRec).^2, 3) ./ ...
                max(TotalWeight, eps('like', TotalWeight));
            uncertaintyRaw = sqrt(residualEnergy ./ residualDenominator);
            uncertaintyRawCPU = gather(real(uncertaintyRaw));
            coreMask = reliableTVMask & confidenceMap >= 0.9;
            if any(coreMask(:))
                uncertaintyBaseline = median(uncertaintyRawCPU(coreMask), 'all');
            else
                uncertaintyBaseline = median(uncertaintyRawCPU(reliableTVMask), 'all');
            end
            uncertaintyExcess = max(uncertaintyRawCPU - uncertaintyBaseline, 0);
            positiveUncertainty = uncertaintyExcess(reliableTVMask & uncertaintyExcess > 0);
            if isempty(positiveUncertainty)
                uncertaintyScale = 1;
            else
                uncertaintyScale = prctile(positiveUncertainty, 95);
            end
            uncertaintyMap = min(uncertaintyExcess ./ max(uncertaintyScale, eps), 1);

            % User-adjustable TV parameters.
            lam_tv_center = lam_tv;
            lam_tv_edge = 50* lam_tv;
            lam_tv_power = 5;
            coverageRiskWeight = 0.8;
            uncertaintyRiskWeight = 0.2;
            tvCenterHeightRatio = 0.65;
            tvCenterWidthRatio = 0.65;
            cropBoundaryRiskPower = 2;
            directionGain = 0.4;

            % Define the protected center rectangle relative to the bounding
            % box of the reliable field, rather than the zero-padded canvas.
            validRows = find(any(reliableTVMask, 2));
            validCols = find(any(reliableTVMask, 1));
            validTop = validRows(1);
            validBottom = validRows(end);
            validLeft = validCols(1);
            validRight = validCols(end);
            validHeight = validBottom - validTop + 1;
            validWidth = validRight - validLeft + 1;

            if isfield(MainPara, 'TV_ROI') && ~isempty(MainPara.TV_ROI)
                validCenterCol = MainPara.TV_ROI(1) + MainPara.TV_ROI(3)/2;
                validCenterRow = MainPara.TV_ROI(2) + MainPara.TV_ROI(4)/2;
                halfCenterWidth = MainPara.TV_ROI(3)/2;
                halfCenterHeight = MainPara.TV_ROI(4)/2;

                [mapColsGrid, mapRowsGrid] = meshgrid(1:size(confidenceMap,2), ...
                    1:size(confidenceMap,1));
                distanceOutsideY = max(abs(mapRowsGrid - validCenterRow) - halfCenterHeight, 0);
                distanceOutsideX = max(abs(mapColsGrid - validCenterCol) - halfCenterWidth, 0);
                availableEdgeY = max(validHeight / 2 - halfCenterHeight, 1);
                availableEdgeX = max(validWidth / 2 - halfCenterWidth, 1);
                normalizedDistanceY = distanceOutsideY ./ availableEdgeY;
                normalizedDistanceX = distanceOutsideX ./ availableEdgeX;
                cropBoundaryRisk = min(max(normalizedDistanceX, normalizedDistanceY), 1) .^ ...
                    cropBoundaryRiskPower;
                centerProtectMask = (distanceOutsideX == 0) & (distanceOutsideY == 0) & ...
                    reliableTVMask;
            else
                cropBoundaryRisk = zeros(size(confidenceMap), 'like', confidenceMap);
                centerProtectMask = false(size(confidenceMap));
            end

            coverageRisk = (1 - confidenceMap).^lam_tv_power;
            uncertaintyRisk = uncertaintyMap .* cropBoundaryRisk;
            riskMap = max(cropBoundaryRisk, coverageRiskWeight .* coverageRisk);
            riskMap = max(riskMap, uncertaintyRiskWeight .* uncertaintyRisk);
            if any(centerProtectMask(:))
                riskMap(centerProtectMask) = 0;
            end
            riskMap(~reliableTVMask) = 0;
            riskMap = min(max(riskMap, 0), 1);
            LambdaMap = lam_tv_center + (lam_tv_edge - lam_tv_center) .* riskMap;
            LambdaMap(~reliableTVMask) = 0;

            % Direction enhancement acts only on the increment above the
            % center lambda, so it cannot increase TV inside the protected ROI.
            sourceM = MNZ_result.M;
            sourceN = MNZ_result.N;
            if isfield(MNZ_result, 'orig_M'), sourceM = MNZ_result.orig_M; end
            if isfield(MNZ_result, 'orig_N'), sourceN = MNZ_result.orig_N; end
            shiftMagnitude = max(abs([sourceM, sourceN]));
            if shiftMagnitude > 0
                directionWeightY = 1 + directionGain * abs(sourceN) / shiftMagnitude;
                directionWeightX = 1 + directionGain * abs(sourceM) / shiftMagnitude;
            else
                directionWeightY = 1;
                directionWeightX = 1;
            end

            validGradVertical = reliableTVMask & reliableTVMask([2:end,end],:);
            validGradVertical(end,:) = false;
            validGradHorizontal = reliableTVMask & reliableTVMask(:,[2:end,end]);
            validGradHorizontal(:,end) = false;
            validGradientMask = cat(3, validGradVertical, validGradHorizontal);
            lambdaBaseVertical = max(LambdaMap, LambdaMap([2:end,end],:));
            lambdaBaseHorizontal = max(LambdaMap, LambdaMap(:,[2:end,end]));
            lambdaVertical = lam_tv_center + directionWeightY .* ...
                max(lambdaBaseVertical - lam_tv_center, 0);
            lambdaHorizontal = lam_tv_center + directionWeightX .* ...
                max(lambdaBaseHorizontal - lam_tv_center, 0);
            LambdaGradient = cat(3, lambdaVertical, lambdaHorizontal);
            LambdaGradient(~validGradientMask) = 0;

            v_est = zeros(size(Object,1), size(Object,2), 2, 'like', Object);
            w_est = zeros(size(Object,1), size(Object,2), 2, 'like', Object);
            LambdaGradientLike = cast(LambdaGradient, 'like', Object);
            validGradientMaskLike = cast(validGradientMask, 'like', Object);
            for subiter = 1:n_subiters_tv
                w_next = v_est + (1 / (8 * gam_tv)) * ...
                    Df(Object - gam_tv * DTf(v_est));
                w_next = min(abs(w_next), LambdaGradientLike) .* exp(1j .* angle(w_next));
                w_next = w_next .* validGradientMaskLike;
                v_est = w_next + (subiter / (subiter + 3)) * (w_next - w_est);
                w_est = w_next;
            end
            Object_tv_candidate = Object - gam_tv * DTf(w_est);
            Object(reliableTVMask) = Object_tv_candidate(reliableTVMask);
            Object(reliableTVMask) = min(abs(Object(reliableTVMask)), 1.05) .* ...
                exp(1j .* angle(Object(reliableTVMask)));

            save(fullfile(foldername, sprintf('AdaptiveTV_iter%04d.mat', iIte)), ...
                'coverageCount', 'confidenceMap', 'uncertaintyRawCPU', ...
                'uncertaintyBaseline', 'uncertaintyMap', 'coverageRisk', ...
                'uncertaintyRisk', 'cropBoundaryRisk', 'centerProtectMask', ...
                'riskMap', 'reliableTVMask', 'LambdaMap', 'LambdaGradient', ...
                'lam_tv_center', 'lam_tv_edge', 'lam_tv_power', ...
                'coverageRiskWeight', 'uncertaintyRiskWeight', ...
                'tvCenterHeightRatio', 'tvCenterWidthRatio', ...
                'cropBoundaryRiskPower', 'directionGain', ...
                'sourceM', 'sourceN', 'directionWeightY', 'directionWeightX');
            imwrite(uint8(confidenceMap * 255), fullfile(foldername, ...
                sprintf('AdaptiveTV_confidence_iter%04d.png', iIte)));
            imwrite(uint8(cropBoundaryRisk * 255), fullfile(foldername, ...
                sprintf('AdaptiveTV_cropBoundaryRisk_iter%04d.png', iIte)));
            imwrite(uint8(centerProtectMask * 255), fullfile(foldername, ...
                sprintf('AdaptiveTV_centerProtectMask_iter%04d.png', iIte)));
            imwrite(uint8(uncertaintyMap * 255), fullfile(foldername, ...
                sprintf('AdaptiveTV_uncertainty_iter%04d.png', iIte)));
            imwrite(uint8(riskMap * 255), fullfile(foldername, ...
                sprintf('AdaptiveTV_risk_iter%04d.png', iIte)));
            imwrite(uint8(LambdaMap ./ max(lam_tv_edge, eps) * 255), ...
                fullfile(foldername, sprintf('AdaptiveTV_lambda_iter%04d.png', iIte)));
        end
        if use_support
            ImgRec = ImgRec .* cast(SupportMask, 'like', ImgRec);
        end

        if isfield(MNZ_result, 'UnionMask')
            outputUnionMask = logical(MNZ_result.UnionMask);
        else
            outputUnionMask = TotalWeight > outputOptions.validMaskThreshold;
        end
        if isfield(MNZ_result, 'TrustedMask')
            outputTrustedMask = logical(MNZ_result.TrustedMask);
        elseif isfield(MNZ_result, 'ValidMaskHard')
            outputCoverageCount = sum(cat(3, MNZ_result.ValidMaskHard{:}), 3);
            outputTrustedMask = outputCoverageCount >= 2;
        else
            outputTrustedMask = outputUnionMask;
        end

        % The main Object remains the maximum union-FOV result.
        [Object, outputUnionMask, outputBounds] = cropToValidFOV( ...
            Object, outputUnionMask, outputOptions);
        if outputOptions.cropToValidFOV
            outputTrustedMask = outputTrustedMask( ...
                outputBounds(2):outputBounds(4), outputBounds(1):outputBounds(3));
        end
        outputValidMask = outputUnionMask; % Backward-compatible variable name.

        % Save an additional stable-area object without changing the main Object.
        Object_trusted = Object;
        Object_trusted(~outputTrustedMask) = 0;
        trustedMinCoverage = 2;
        save(fullfile(foldername, 'OutputFOV.mat'), ...
            'outputValidMask', 'outputUnionMask', 'outputTrustedMask', ...
            'outputBounds', 'trustedMinCoverage');
        save(fullfile(foldername, sprintf("Object_iter%04d.mat", iIte)), 'Object');
        save(fullfile(foldername, sprintf("Object_trusted_iter%04d.mat", iIte)), ...
            'Object_trusted', 'outputTrustedMask', 'trustedMinCoverage');
        % <<< OUTPUT-UNION-TRUSTED-FOV-END (Codex)


        % 3. 外围盲区相位置零
        % Invalid pixels were zero-filled by cropToValidFOV above.
        % ---------------------------------------------

        % Save amplitude / phase images
        ampObj = abs(Object);
        phsObj = angle(Object);


        % 使用 imadjust 拉伸对比度
        % >>> AMP-DISPLAY-FIX-BEGIN (Codex)
        amp_gray = mat2gray(ampObj);
        % <<< AMP-DISPLAY-FIX-END (Codex)
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


        imwrite(amp_gray_stretched, fullfile(foldername, sprintf("Object_amp_iter%04d.png", save_iter)));
        imwrite(phs_gray_stretched, fullfile(foldername, sprintf("Object_phs_iter%04d.png", save_iter)));

        count = count + 1;
    end

    % --- 如果判定已收敛，跳出迭代循环 ---
    if exist('is_converged', 'var') && is_converged
        break;
    end
end


end

function [field, validMask, bounds] = cropToValidFOV(field, validMask, options)
% Remove replicated padding without discarding synthetic-aperture content.
if ~any(validMask(:))
    error('Output valid mask is empty; check illumination shifts and padding.');
end

if options.cropToValidFOV
    [rows, cols] = find(validMask);
    bounds = [min(cols), min(rows), max(cols), max(rows)]; % [x1 y1 x2 y2]
    field = field(bounds(2):bounds(4), bounds(1):bounds(3));
    validMask = validMask(bounds(2):bounds(4), bounds(1):bounds(3));
else
    bounds = [1, 1, size(field, 2), size(field, 1)];
end

if options.zeroFillInvalid
    field(~validMask) = 0;
end
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
