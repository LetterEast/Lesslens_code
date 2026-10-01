function outputField = propagateAngularSpectrum(inputField, pixelSize, wavelength, distance)
%PROPAGATE Band-limited angular spectrum propagation.
% pixelSize [m/pixel], wavelength [m], distance [m]. Positive distance
% propagates forward; negative distance back propagates toward the sample.
% Uses an available GPU and otherwise computes on the CPU.

if distance == 0
    outputField = inputField;
    return;
end

[numRows, numCols] = size(inputField);

useGPU = exist('gpuDeviceCount', 'file') == 2 && ...
    license('test', 'Distrib_Computing_Toolbox') && gpuDeviceCount > 0;

if useGPU && ~isa(inputField, 'gpuArray')
    inputField = gpuArray(inputField);
end

%% ---------------------------------------------------------
% Spatial-frequency coordinates
% ----------------------------------------------------------
fx = ifftshift( ...
    (-floor(numCols/2):ceil(numCols/2)-1) ...
    / (numCols * pixelSize));

fy = ifftshift( ...
    (-floor(numRows/2):ceil(numRows/2)-1).' ...
    / (numRows * pixelSize));

if useGPU
    fx = gpuArray(fx);
    fy = gpuArray(fy);
end

[FX, FY] = meshgrid(fx, fy);

%% ---------------------------------------------------------
% Angular-spectrum transfer function
% ----------------------------------------------------------
rootArgument = ...
    1 ...
    - (wavelength * FX).^2 ...
    - (wavelength * FY).^2;

propagatingMask = rootArgument >= 0;

H = zeros(size(rootArgument), 'like', inputField);

H(propagatingMask) = exp( ...
    1i * 2*pi/wavelength * distance .* ...
    sqrt(rootArgument(propagatingMask)) );

%% ---------------------------------------------------------
% Band-limited ASM
% ----------------------------------------------------------

Lx = numCols * pixelSize;
Ly = numRows * pixelSize;

fxLimit = 1 ./ ...
    (wavelength * sqrt(1 + (2*abs(distance)/Lx).^2));

fyLimit = 1 ./ ...
    (wavelength * sqrt(1 + (2*abs(distance)/Ly).^2));

% Nyquist limits
fxNyquist = 1 / (2 * pixelSize);
fyNyquist = 1 / (2 * pixelSize);

fxLimit = min(fxLimit, fxNyquist);
fyLimit = min(fyLimit, fyNyquist);

bandMask = ...
    abs(FX) <= fxLimit & ...
    abs(FY) <= fyLimit;

H = H .* bandMask;

%% ---------------------------------------------------------
% Propagation
% ----------------------------------------------------------
outputField = ifft2(fft2(inputField) .* H);

if useGPU
    wait(gpuDevice);
    outputField = gatherIfNeeded(outputField);
end

end
