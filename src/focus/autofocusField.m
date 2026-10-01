function [bestDistance, diagnostics] = autofocusField(field, pixelSize, wavelength, settings)
%AUTOFOCUSFIELD Scan sample distance [m] with a selectable sharpness metric.
% method='fft' preserves the original full-field high-pass amplitude metric.
% method='adfrft' evaluates intensity in a fractional Fourier domain. Both
% use the same distance grid and the existing reconstructed complex field.
validateattributes(settings.prior, {'numeric'}, {'real','finite','scalar'});
validateattributes(settings.halfRange, {'numeric'}, {'real','finite','scalar','nonnegative'});
validateattributes(settings.step, {'numeric'}, {'real','finite','scalar','positive'});
method = 'fft';
if isfield(settings,'method'), method = validatestring(settings.method, {'fft','adfrft'}); end
frftOrder = [];
cropSize = 0;
metricName = 'High-pass amplitude FFT magnitude sum';
methodLabel = 'FFT high-pass amplitude';
if strcmp(method,'adfrft')
    frftOrder = 0.8;
    if isfield(settings,'frftOrder'), frftOrder = settings.frftOrder; end
    validateattributes(frftOrder, {'numeric'}, {'real','finite','scalar','positive','<=',1});
    cropSize = 256;
    if isfield(settings,'cropSize'), cropSize = settings.cropSize; end
    validateattributes(cropSize, {'numeric'}, {'real','finite','scalar','integer','nonnegative'});
    metricName = 'mean(abs(FrFT_p(abs(U).^2))) + (p-1)*mean(abs(U).^2)';
    methodLabel = sprintf('ADFrFT (p = %.3g)',frftOrder);
end
% Crop after full-field propagation: do not introduce a new aperture into
% propagation or discard the phase recovered by APRW.
[height,width] = size(field);
roiBounds = [1,1,width,height];
roiMode = 'full';
if strcmp(method,'adfrft')
    roiMode = 'center';
    if isfield(settings,'roiMode')
        roiMode = validatestring(settings.roiMode,{'center','full','manual'});
    end
end
if strcmp(roiMode,'manual')
    if isfield(settings,'roiBounds') && ~isempty(settings.roiBounds)
        roiBounds = settings.roiBounds;
    else
        preview = propagateAngularSpectrum(field,pixelSize,wavelength,-settings.prior);
        roiBounds = selectFocusROI(abs(preview));
    end
    validateattributes(roiBounds,{'numeric'}, ...
        {'real','finite','integer','positive','vector','numel',4});
    roiBounds = double(roiBounds(:).');
    if roiBounds(3)>width || roiBounds(4)>height || ...
            roiBounds(3)<=roiBounds(1) || roiBounds(4)<=roiBounds(2)
        error('Lensless:InvalidFocusROI','ROI must lie within the field and be at least 2-by-2 pixels.');
    end
elseif strcmp(roiMode,'center') && cropSize > 0
    count = min([cropSize,height,width]);
    firstRow = floor((height-count)/2)+1;
    firstColumn = floor((width-count)/2)+1;
    roiBounds = [firstColumn,firstRow,firstColumn+count-1,firstRow+count-1];
else
    roiMode = 'full';
end
distances = settings.prior + (-settings.halfRange:settings.step:settings.halfRange);
if isempty(distances), distances = settings.prior; end
focusMetric = zeros(size(distances));
for index = 1:numel(distances)
    object = propagateAngularSpectrum(field, pixelSize, wavelength, -distances(index));
    if strcmp(method,'adfrft')
        roi = object(roiBounds(2):roiBounds(4),roiBounds(1):roiBounds(3));
        intensity = double(gatherIfNeeded(abs(roi).^2));
        transformed = fractionalFourier2(intensity,frftOrder);
        focusMetric(index) = mean(abs(transformed),'all') + ...
            (frftOrder-1)*mean(intensity,'all');
    else
        amplitude = abs(object);
        highPass = amplitude - imgaussfilt(amplitude, 2);
        focusMetric(index) = double(gatherIfNeeded(sum(abs(fft2(highPass)), 'all')));
    end
end
if any(~isfinite(focusMetric))
    error('Lensless:InvalidFocusMetric', 'Autofocus produced non-finite sharpness scores.');
end
[bestMetric, bestIndex] = max(focusMetric);
bestDistance = distances(bestIndex);
diagnostics = struct('distances', distances, 'focusMetric', focusMetric, ...
    'bestDistance', bestDistance, 'bestIndex', bestIndex, ...
    'bestMetric', bestMetric, 'prior', settings.prior, ...
    'halfRange', settings.halfRange, 'step', settings.step, ...
    'method', method, 'methodLabel', methodLabel, 'metricName', metricName, ...
    'frftOrder', frftOrder, 'cropSize', cropSize, 'roiBounds', roiBounds, 'roiMode', roiMode, ...
    'searchPerformed', numel(distances) > 1, ...
    'peakAtBoundary', numel(distances) > 1 && ...
        (bestIndex == 1 || bestIndex == numel(distances)), ...
    'flatCurve', numel(distances) > 1 && ...
        max(focusMetric)-min(focusMetric) <= eps(max(1, max(abs(focusMetric))))*16);
end
