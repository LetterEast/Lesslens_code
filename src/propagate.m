function outputField = propagate(inputField, pixelSize, wavelength, distance)
%PROPAGATE On-axis angular-spectrum propagation with optional GPU acceleration.
%   Positive distance propagates forward; negative distance propagates back.

if distance == 0
    outputField = inputField;
    return;
end

[numRows, numCols] = size(inputField);
useGPU = gpuDeviceCount > 0;
if useGPU && ~isa(inputField, 'gpuArray')
    inputField = gpuArray(inputField);
end

frequencyX = ifftshift((-floor(numCols/2):ceil(numCols/2)-1) / ...
    (numCols * pixelSize));
frequencyY = ifftshift((-floor(numRows/2):ceil(numRows/2)-1).' / ...
    (numRows * pixelSize));
if useGPU
    frequencyX = gpuArray(frequencyX);
    frequencyY = gpuArray(frequencyY);
end
[frequencyX, frequencyY] = meshgrid(frequencyX, frequencyY);

rootArgument = 1 - (wavelength * frequencyX).^2 - ...
    (wavelength * frequencyY).^2;
propagating = rootArgument >= 0;
transfer = zeros(size(rootArgument), 'like', inputField);
transfer(propagating) = exp(1i * (2*pi/wavelength) * distance .* ...
    sqrt(rootArgument(propagating)));

outputField = ifft2(fft2(inputField) .* transfer);
if useGPU
    wait(gpuDevice);
    outputField = gather(outputField);
end
end
