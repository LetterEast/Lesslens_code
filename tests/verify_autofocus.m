function verify_autofocus()
%VERIFY_AUTOFOCUS Check degenerate scans and validation without real datasets.
projectRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(projectRoot,'src'));
settings = struct('prior',1e-3,'halfRange',0,'step',1e-5);
[distance, info] = autofocusField(ones(16,16), 3e-6, 514e-9, settings);
assert(distance==settings.prior && ~info.searchPerformed && ~info.peakAtBoundary);
settings.halfRange = 2e-5;
[distance, info] = autofocusField(ones(16,16), 3e-6, 514e-9, settings);
assert(numel(info.distances)==5 && info.flatCurve);
assert(distance==info.distances(1)); % Preserve the original first-maximum rule.
settings.step = 0;
rejected = false;
try
    autofocusField(ones(16,16), 3e-6, 514e-9, settings);
catch
    rejected = true;
end
assert(rejected, 'A zero scan step must fail before propagation.');
% ADFrFT on a fixed distance must match an independently evaluated formula.
settings = struct('prior',0,'halfRange',0,'step',1e-5, ...
    'method','adfrft','frftOrder',0.8,'cropSize',12);
[x,y] = meshgrid(1:22,1:18);
field = (0.5+0.1*cos(x/3).*sin(y/2)).*exp(1i*x/7);
[distance,info] = autofocusField(field,3e-6,514e-9,settings);
intensity = abs(field(4:15,6:17)).^2;
expected = mean(abs(fractionalFourier2(intensity,0.8)),'all')-0.2*mean(intensity,'all');
assert(distance==0 && strcmp(info.method,'adfrft') && info.frftOrder==0.8);
assert(isequal(info.roiBounds,[6,4,17,15]));
assert(abs(info.focusMetric-expected)<1e-12);
settings.halfRange = 2e-5;
[~,info] = autofocusField(field,3e-6,514e-9,settings);
assert(info.searchPerformed && numel(info.distances)==5 && all(isfinite(info.focusMetric)));
assert(info.bestDistance==info.distances(info.bestIndex));
settings.cropSize = 0;
[~,info] = autofocusField(field,3e-6,514e-9,settings);
assert(isequal(info.roiBounds,[1,1,22,18]));
settings.method = 'invalid';
rejected = false;
try
    autofocusField(field,3e-6,514e-9,settings);
catch
    rejected = true;
end
assert(rejected, 'An unknown method must not silently select another metric.');
settings.method = 'adfrft';
settings.prior = 0;
settings.halfRange = 0;
settings.roiMode = 'manual';
settings.roiBounds = [2,3,10,9];
[~,info] = autofocusField(field,3e-6,514e-9,settings);
intensity = abs(field(3:9,2:10)).^2;
expected = mean(abs(fractionalFourier2(intensity,0.8)),'all')-0.2*mean(intensity,'all');
assert(isequal(info.roiBounds,settings.roiBounds) && strcmp(info.roiMode,'manual'));
assert(abs(info.focusMetric-expected)<1e-12);
settings.roiBounds = [2,3,100,9];
rejected = false;
try
    autofocusField(field,3e-6,514e-9,settings);
catch exception
    rejected = strcmp(exception.identifier,'Lensless:InvalidFocusROI');
end
assert(rejected,'Out-of-image manual bounds must be rejected.');
fprintf('Autofocus edge-case checks passed.\n');
end
