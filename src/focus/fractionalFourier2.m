function output = fractionalFourier2(input, order)
%FRACTIONALFOURIER2 Separable 2-D fractional Fourier transform.
% Chirp/interpolation/convolution algorithm of Ozaktas et al., IEEE TSP 44,
% 2141 (1996). Uses sinc interpolation and a centered sampling grid.
% Order 0 is identity; order 1 is a unitary centered
% DFT. Noninteger orders are numerical approximations on a finite grid.
validateattributes(input, {'single','double'}, {'2d','nonempty','finite'});
validateattributes(order, {'numeric'}, {'real','scalar','finite'});
assert(all(size(input)>=2), 'FrFT requires at least two samples on each axis.');
output = transformRows(double(gatherIfNeeded(input)), order);
output = transformRows(output.', order).';
end

function output = transformRows(input, order)
n = size(input,2);
order = mod(order,4);
if order == 0, output = input; return; end
if order == 2, output = fliplr(input); return; end
if order == 1, output = centeredDFT(input,false); return; end
if order == 3, output = centeredDFT(input,true); return; end
if order > 2
    order = order-2;
    input = fliplr(input);
end
if order > 1.5
    order = order-1;
    input = centeredDFT(input,false);
end
if order < 0.5
    order = order+1;
    input = centeredDFT(input,true);
end
% Band-limited 2x interpolation, without a Signal Processing Toolbox dependency.
upsampled = complex(zeros(size(input,1),2*n-1));
upsampled(:,1:2:end) = input;
t = (-(2*n-3):(2*n-3))/2;
kernel = ones(size(t));
nonzero = t~=0;
kernel(nonzero) = sin(pi*t(nonzero))./(pi*t(nonzero));
interpolated = linearConvolution(upsampled,kernel);
interpolated = interpolated(:,2*n-2:4*n-4);
alpha = order*pi/2;
k = -n:n-2;
chirp = exp(-1i*pi/n*tan(alpha/2)/4*k.^2);
c = pi/n/sin(alpha)/4;
lags = -(2*n-2):2*n-2;
convolved = linearConvolution(chirp.*interpolated,exp(1i*c*lags.^2));
output = convolved(:,2*n-1:4*n-3)*sqrt(c/pi);
output = chirp.*output;
output = exp(-1i*(1-order)*pi/4)*output(:,1:2:end);
end

function output = linearConvolution(input,kernel)
count = size(input,2)+numel(kernel)-1;
fftLength = 2^nextpow2(count);
output = ifft(fft(input,fftLength,2).*fft(kernel,fftLength,2),[],2);
output = output(:,1:count);
end

function output = centeredDFT(input,inverse)
n = size(input,2);
indices = mod((0:n-1)+floor(n/2),n)+1;
output = complex(zeros(size(input)));
if inverse
    output(:,indices) = ifft(input(:,indices),[],2)*sqrt(n);
else
    output(:,indices) = fft(input(:,indices),[],2)/sqrt(n);
end
end
