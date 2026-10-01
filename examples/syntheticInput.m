function inputData = syntheticInput(wavelength)
%SYNTHETICINPUT Small self-contained fixture for examples and regression tests.
% This is a demonstration of the prepared input format, not a measured dataset
% or a validation of reconstruction quality on a physical experiment.
if nargin < 1, wavelength = 514e-9; end
imageSize = [48, 64];
g = struct('Z',0.05,'M',0,'N',0,'orig_M',0,'orig_N',0, ...
    'padSize',[8,8],'orig_size',[32,48]);
mask = false(imageSize);
mask(9:40,9:56) = true;
g.ValidMaskHard = {mask, mask};
[x,y] = meshgrid(1:imageSize(2), 1:imageSize(1));
intensity = 0.4 + 0.15*sin(x/3).*cos(y/4) + 0.05*cos((x+y)/2);
inputData = struct('images', {{intensity, intensity*0.9}}, ...
    'geometry', g, 'distanceSteps', [0,1e-3], ...
    'wavelength', wavelength, 'pixelSize', 3e-6);
end
