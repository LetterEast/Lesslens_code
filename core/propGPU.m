function newU = propGPU(U, sampling_rate, lambda, z, x0, y0, validSize)
% =========================================================================
% Copyright (c) 2026 [JiNan University]
%
% Implementation of the Shifted Angular Spectrum Method (Shift-AS).
% This code is a numerical reproduction based on the formulation by:
%
% Kyoji Matsushima, "Shifted angular spectrum method for off-axis numerical
% propagation," Optics Express, Vol. 18, Issue 17, pp. 18453-18463 (2010).
%
% Key features implemented:
% 1. Off-axis coordinate shift using Fourier Shift Theorem (Eq. 7 & 8)[cite: 147, 149].
% 2. Band-limited mask to avoid aliasing in shifted propagation[cite: 192].
% 3. Frequency-domain demodulation/remodulation for spectrum centering.
%
% Licensed under the MIT License.
% =========================================================================

[M, N] = size(U);
if nargin < 5 || isempty(x0), x0 = 0; end
if nargin < 6 || isempty(y0), y0 = 0; end
if nargin < 7 || isempty(validSize), validSize = [M, N]; end
x0 = 0;
validateattributes(validSize, {'numeric'}, ...
    {'vector', 'numel', 2, 'positive', 'finite'});
validSize = min(round(validSize(:).'), [M, N]);
useGPU = (gpuDeviceCount > 0);
if useGPU && ~isa(U, 'gpuArray')
    U = gpuArray(U);
end

% --- 1. 生成空间中心坐标网格 ---
x_vec = ((1:N) - floor(N/2) - 1) * sampling_rate;
y_vec = ((1:M) - floor(M/2) - 1)' * sampling_rate;
if useGPU
    x_vec = gpuArray(x_vec);
    y_vec = gpuArray(y_vec);
end
[X, Y] = meshgrid(x_vec, y_vec);

% --- 2. 解析计算理论载波频率 ---
L = sqrt(z^2 + x0^2 + y0^2);

if L == 0
    fx = 0;
    fy = 0;
else
    fx = (x0 / L) / lambda * sign(z);
    fy = (y0 / L) / lambda * sign(z);
end
% --- 3. 完美解调 (Demodulation) ---
U_base = U .* exp(-1i * 2 * pi * (fx * X + fy * Y));

% --- 4. 空间域边缘平滑 ---
win_width = 0;
W_spatial = tukeywin(M, win_width) * tukeywin(N, win_width)';
if useGPU, W_spatial = gpuArray(W_spatial); end
U_base = U_base .* W_spatial;

% --- 5. 边缘复制扩充 (替代原本的补零) ---
% 纯补零会在边缘产生巨大的阶跃（从背景亮度突然掉到 0），传播时会产生强烈的菲涅尔衍射波纹（吉布斯振铃）侵入原图边缘
% 边缘复制扩充（Replicate Padding）可以完美避免这种阶跃
pad_factor = 1; % 增加补零系数，防止大衍射角下的卷边问题
pad_M = M * pad_factor;
pad_N = N * pad_factor;

start_row = floor((pad_M - M)/2) + 1;
start_col = floor((pad_N - N)/2) + 1;

% --- 方案 A: 纯零填充 (测试用) ---
% U_pad = zeros(pad_M, pad_N, 'like', U_base);
% U_pad(start_row : start_row+M-1, start_col : start_col+N-1) = U_base;

% --- 方案 B: 边缘复制扩充 (原本的代码，效果好可随时切回) ---
idx_r = [ones(1, start_row - 1), 1:M, M * ones(1, pad_M - start_row - M + 1)];
idx_c = [ones(1, start_col - 1), 1:N, N * ones(1, pad_N - start_col - N + 1)];
U_pad = U_base(idx_r, idx_c);

% --- 6. 频率坐标网格 ---
du = 1 / (pad_N * sampling_rate);
dv = 1 / (pad_M * sampling_rate);
u_vec = ((1:pad_N) - floor(pad_N/2) - 1) * du;
v_vec = ((1:pad_M) - floor(pad_M/2) - 1) * dv;
if useGPU
    u_vec = gpuArray(u_vec);
    v_vec = gpuArray(v_vec);
end
[U_freq, V_freq] = meshgrid(u_vec, v_vec);

% --- 7. 前向 FFT (此时频谱绝对完美居中) ---
Uspectrum = fftshift(fft2(ifftshift(U_pad)));

% 计算真实物理频率 (Demodulation limits spectrum to 0, so true frequency is shifted back)
U_true = U_freq + fx;
V_true = V_freq + fy;

% --- 8. 标准轴向角谱传递函数 ---
w_sq = 1 - (lambda * U_true).^2 - (lambda * V_true).^2;
valid_w = w_sq > 0;
H_AS = zeros(size(w_sq), 'like', w_sq);
H_AS(valid_w) = exp(1i * (2 * pi / lambda) * z * sqrt(w_sq(valid_w)));

% 软切断带限掩模 (Matsushima 2010 严格带限 Shifted-AS mask)
% The band limit is determined by the measured/valid FOV, not by the
% artificial computational padding around it.
X_len = validSize(2) * sampling_rate;
Y_len = validSize(1) * sampling_rate;

u1 = (-x0 + X_len/2) / lambda / sqrt(z^2 + (-x0 + X_len/2)^2);
u2 = (-x0 - X_len/2) / lambda / sqrt(z^2 + (-x0 - X_len/2)^2);
f_nyquist = 1 / (2 * sampling_rate);
u_max = min(max(u1, u2), f_nyquist);
u_min = max(min(u1, u2), -f_nyquist);

v1 = (-y0 + Y_len/2) / lambda / sqrt(z^2 + (-y0 + Y_len/2)^2);
v2 = (-y0 - Y_len/2) / lambda / sqrt(z^2 + (-y0 - Y_len/2)^2);
v_max = min(max(v1, v2), f_nyquist);
v_min = max(min(v1, v2), -f_nyquist);

if abs(z) <= eps(class(z))
    mask = ones(size(U_true), 'like', U_true);
else
    mask = softRectBandlimit( ...
        U_true, V_true, u_min, u_max, v_min, v_max, 0.1);
end
mask = ones(size(U_true), 'like', U_true);
% --- 9. 傅里叶平移定理 (核心物理转移) ---
H_shift = exp(-1i * 2 * pi * (U_true * x0 + V_true * y0));

% --- 10. 频域滤波与逆变换 ---
newUspectrum = Uspectrum .* H_AS .* H_shift .* mask;
newU_pad = fftshift(ifft2(ifftshift(newUspectrum)));

% --- 11. 截取原始物理视场 ---
newU_base = newU_pad(start_row : start_row+M-1, start_col : start_col+N-1);

% --- 12. 重新调制 (Remodulation) ---
newU = newU_base .* exp(1i * 2 * pi * (fx * X + fy * Y));

if useGPU
    wait(gpuDevice);
    newU = gather(newU);
end
end
function mask = softRectBandlimit(U, V, uMin, uMax, vMin, vMax, rollFraction)
% Separable Matsushima passband with a cosine transition at each boundary.
maskU = softBand1D(U, uMin, uMax, rollFraction);
maskV = softBand1D(V, vMin, vMax, rollFraction);
mask = maskU .* maskV;
end

function weight = softBand1D(freq, fMin, fMax, rollFraction)
weight = zeros(size(freq), 'like', freq);
bandwidth = fMax - fMin;
if bandwidth <= 0
    return;
end
rollWidth = max(bandwidth * rollFraction, eps(class(bandwidth)));
innerMin = min(fMin + rollWidth, (fMin + fMax) / 2);
innerMax = max(fMax - rollWidth, (fMin + fMax) / 2);
weight(freq >= innerMin & freq <= innerMax) = 1;
idx = freq >= fMin & freq < innerMin;
weight(idx) = 0.5 * (1 - cos(pi * (freq(idx) - fMin) / ...
    max(innerMin - fMin, eps(class(bandwidth)))));
idx = freq > innerMax & freq <= fMax;
weight(idx) = 0.5 * (1 + cos(pi * (freq(idx) - innerMax) / ...
    max(fMax - innerMax, eps(class(bandwidth)))));
end
