function verify_frft()
%VERIFY_FRFT Exact integer orders and approximate Gaussian invariance.
projectRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(projectRoot,'src'));
[x,y] = meshgrid(1:20,1:16);
input = sin(x/3).*cos(y/4)+1i*cos((x+y)/5);
assert(isequal(fractionalFourier2(input,0),input));
expected = fftshift(fft2(ifftshift(input)))/sqrt(numel(input));
actual = fractionalFourier2(input,1);
assert(norm(actual-expected,'fro')/norm(expected,'fro')<1e-12);
assert(isequal(fractionalFourier2(input,2),flipud(fliplr(input))));
assert(norm(fractionalFourier2(actual,3)-input,'fro')/norm(input,'fro')<1e-12);
[x,y] = meshgrid(((0:127)-64)/sqrt(128),((0:95)-48)/sqrt(96));
gaussian = exp(-pi*(x.^2+y.^2));
for order = [0.2,0.8,1.3]
    transformed = fractionalFourier2(gaussian,order);
    assert(norm(abs(transformed)-gaussian,'fro')/norm(gaussian,'fro')<0.02);
end
fprintf('Fractional Fourier transform checks passed.\n');
end
