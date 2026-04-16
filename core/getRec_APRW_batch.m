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
%   predefined_rect    : predefined cropping region for calibration, to bypass interactive selection.
%
% Output:
%   Rec_nInterative    : cell, recorded complex field ImgRec at each record point
%   rect               : the bounding box used for cropping

%% ---------- Basic checks & record setting ----------
nIterative = MainPara.nIterative;

if nargin < 4 || isempty(iIte_record)
    iIte_record = nIterative; % default: only record last iteration
end

if nargin < 5
    predefined_rect = [];
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

%% ---------- Weighted feedback coefficients ----------

a = 0.7;
b = -0.13*a + 0.6;   %

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

%% ---------- Precalculate alignment shifts (Data-Driven HSSA) ----------
% Align images between different Z heights with the SAME illumination angle,
% and set the original M/N geometric shift calculation to 0.
num_xy = sum(zSet == 0);
if num_xy == 0
    num_xy = 1;
end
sa_shift = zeros(n_img, 2); % [row_shift, col_shift]
for iImg = 1:n_img
    ref_idx = mod(iImg - 1, num_xy) + 1; % Corresponding image in Z1 block
    if iImg ~= ref_idx
        % Measure subpixel shift relative to its corresponding Z1 reference
        [out, ~] = dftregistration(fft2(img_set{ref_idx}), fft2(img_set{iImg}), 100);
        sa_shift(iImg, :) = [out(3), out(4)];
    end
end

%% ---------- Create a unique result root folder ----------
runStamp = datestr(now, 'yyyymmdd_HHMMSS');
Parent_foldername = fullfile(".\ResultFolder", sprintf("APRW_%s_N%d_I%d_rec%d", runStamp, n_img, nIterative, iIte_record));
if ~exist(Parent_foldername, 'dir'); mkdir(Parent_foldername); end

%% ---------- Initialization ----------
% Initialize complex field estimate at plane 1 using sqrt of first intensity
ImgRec = sqrt(img_set{1});
ImgRec = ImgRec .* IllumPattern;  % apply illumination modulation if needed

% Store modulated guesses from previous iterations for weighted feedback
ImgRec_record1 = cell(1, (n_img-1));  % g~_{k-1}(n)
ImgRec_record2 = cell(1, (n_img-1));  % g~_{k-2}(n)

%% ---------- Visualization ----------
hFig = figure('Name', 'APRW / Weighted Feedback Monitor');
ax1 = subplot(1,2,1, 'Parent', hFig);
ax2 = subplot(1,2,2, 'Parent', hFig);

% placeholders
hAmp = imshow(zeros(mRow,nCol), [], 'Parent', ax1); title(ax1, 'Amplitude');
hPhs = imshow(zeros(mRow,nCol), [], 'Parent', ax2); title(ax2, 'Phase');

%% ---------- Reconstruction process ----------
tic;
DistanceInterval = 0;
zSet = zeros(n_img,1);
for iIte = 1 : nIterative

    % One iteration: build object guesses for planes 2..N
    ImgRec_record = cell(1, (n_img-1));

    for iImg = 2 : n_img
        if iIte == 1
            Interval = DistanceIntervalSet(iImg);
            DistanceInterval = Interval + DistanceInterval;
            zSet(iImg) = DistanceInterval;
        end


        % Current propagation distance from plane#1 to plane#iImg
        z = zSet(iImg);      % [m]
        z_mm = z * 1e3;      % [mm] for display

        fprintf('Plane 1 -> %d propagation distance: %0.3f mm\n', iImg, z_mm);

        % (1) Forward propagate to detector plane iImg
        x = -MNZ_result.M*PixelSize*(z/D_LED2CCD);
        y = -MNZ_result.N*PixelSize*(z/D_LED2CCD);
        lightOnDetector = propGPU(ImgRec, PixelSize, WaveLength, z, x, y);

        % (2) Amplitude replacement: measured amplitude = sqrt(intensity)
        thisAmpImg      = img_set{iImg}.^0.5;
        lightOnDetector = thisAmpImg .* exp(1j .* angle(lightOnDetector));

        % (3) Back propagate to object plane => g_k(n)
        ImgRec_record{iImg-1} = propGPU(lightOnDetector, PixelSize, WaveLength, -z, -x, -y);

        % (4) Weighted feedback (start when k > 2)
        if iIte <= 2
            g = ImgRec_record{iImg-1}; % no feedback for first two iterations
        else
            g = (1 + a + b) .* ImgRec_record{iImg-1} ...
                - a .* ImgRec_record1{iImg-1} ...
                - b .* ImgRec_record2{iImg-1};
        end
        ImgRec_record{iImg-1} = g;

        % ----- Display amplitude & phase for current plane guess -----
        roi = g(MarginY+1:MarginY+mRow, MarginX+1:MarginX+nCol);

        AmpImgRec = gather(abs(roi));
        PhsImgRec = gather(angle(roi));

        AmpImgRec = mat2gray(AmpImgRec);
        PhsImgRec = mat2gray(PhsImgRec);
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

    Stack = cat(3, ImgRec_record1{:});
    ImgRec = mean(Stack, 3);

    %% ---------- Record & save results ----------
    if (mod(iIte, iIte_record) == 0) && (count <= record_num)

        % Subfolder for this recorded iteration
        foldername = fullfile(Parent_foldername, sprintf("Iter_%04d", iIte));
        if ~exist(foldername, 'dir'); mkdir(foldername); end

        % Store complex field (at reconstruction plane) in output cell
        Rec_nInterative{count} = ImgRec;

        % calibrated sample-to-CCD distance
        % D_Sample2CCD = 0.00184; % [m]

        if iIte == iIte_record
            if isempty(predefined_rect)
                Ishow = gather(mat2gray(abs(ImgRec)));
                figure, [~, rect] = imcrop(Ishow); close
            else
                rect = predefined_rect;
            end
        end
        rect = round(rect);                      % [x y w h]
        x1 = max(1, floor(rect(1)) + 1);
        y1 = max(1, floor(rect(2)) + 1);
        w  = max(1, rect(3));
        h  = max(1, rect(4));

        x2 = min(size(ImgRec,2), x1 + w - 1);
        y2 = min(size(ImgRec,1), y1 + h - 1);

        ImgRec_crop = ImgRec(y1:y2, x1:x2);
        D_Sample2CCD = D_Sample2CCD_Calibration(ImgRec_crop, ParaD_Sample2CCD, MainPara, foldername);
        %
        % D_Sample2CCD = D_Sample2CCD_Calibration(ImgRec, ParaD_Sample2CCD, MainPara, foldername);

        % Save metadata for reproducibility
        metaPath = fullfile(foldername, "meta.mat");
        save(metaPath, 'MainPara', 'DistanceIntervalSet', 'zSet', 'a', 'b', 'D_Sample2CCD', 'iIte','ImgRec');

        % Back to object plane and save Object
        Object = back2Object(ImgRec, D_Sample2CCD, MainPara);
        save(fullfile(foldername, sprintf("Object_iter%04d.mat", iIte)), 'Object');

        % Save amplitude / phase images
        ampObj = abs(Object);
        phsObj = angle(Object);

        imwrite(imadjust(mat2gray(ampObj)), fullfile(foldername, sprintf("Object_amp_iter%04d.png", iIte)));
        imwrite(mat2gray(phsObj),          fullfile(foldername, sprintf("Object_phs_iter%04d.png", iIte)));

        % Quick view
        hFigObj = figure('Name', sprintf('Object @ Iter %d', iIte));
        subplot(1,2,1), imshow(ampObj, []), xlabel("Object_amp"), axis image
        subplot(1,2,2), imshow(mat2gray(phsObj)), xlabel("Object_phs"), axis image
        drawnow;
        close(hFigObj);

        count = count + 1;
    end
end

toc;
disp('finish reconstruction');
end
