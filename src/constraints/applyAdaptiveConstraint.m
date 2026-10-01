function [field, state, masks] = applyAdaptiveConstraint(field, inputData, settings, state)
%APPLYADAPTIVECONSTRAINT Relax absorption/phase supports in the sample plane.
% Independent implementation of the morphology-based update described by
% Xu et al., Scientific Reports 13, 10267 (2023), doi:10.1038/s41598-023-37423-3.
% Adds relaxation and spherical-illumination normalization to that method.
% settings.distance must be resolved by the caller (metres). whiteAmplitude
% is the incident amplitude at the reference detector, in input data units.
% state belongs to ONE reconstruction run; use [] when input/settings change.
% Zero strength is an exact bypass, including propagation and normalization.
if nargin < 4, state = []; end
masks = struct();
validateattributes(settings.strength, {'numeric'}, ...
    {'real','finite','scalar','>=',0,'<=',1});
if settings.strength == 0, return; end
if isempty(state)
    g = inputData.geometry;
    z = settings.distance;
    validateattributes(z, {'numeric'}, {'real','finite','scalar','<',g.Z});
    white = settings.whiteAmplitude;
    if isempty(white)
        values = inputData.images{1};
        if isfield(g,'ValidMaskHard'), values = values(logical(g.ValidMaskHard{1})); end
        white = sqrt(prctile(double(gatherIfNeeded(max(values(:),0))),99.9));
    end
    validateattributes(white, {'numeric'}, {'real','finite','scalar','positive'});
    support = sampleGeometrySupport(size(field),g,cumsum(inputData.distanceSteps),z);
    [r,c] = find(support.mask);
    state.rows = min(r):max(r); state.cols = min(c):max(c);
    state.valid = support.mask(state.rows,state.cols);
    edgeWidth = 32;
    if isfield(settings,'edgeWidth'), edgeWidth = settings.edgeWidth; end
    validateattributes(edgeWidth, {'numeric'}, {'real','finite','scalar','positive'});
    % Include a false border so a footprint touching the array boundary also
    % tapers. Only the UPDATE is tapered; the reconstructed field is not windowed.
    distance = bwdist(~padarray(state.valid,[1,1],false,'both'));
    distance = max(double(distance(2:end-1,2:end-1))-1,0);
    taper = sin((pi/2)*min(distance/edgeWidth,1)).^2;
    % Smooth coverage before weighting: integer plane-count steps must not
    % become new seams. Low-coverage regions still receive at least half
    % strength away from the exterior boundary.
    coverage = double(support.coverage(state.rows,state.cols));
    sigma = max(edgeWidth/2,1);
    smoothCoverage = imgaussfilt(coverage,sigma,'Padding',0)./ ...
        max(imgaussfilt(double(state.valid),sigma,'Padding',0),eps);
    confidence = min(max(smoothCoverage/max(coverage(:)),0),1);
    state.updateWeight = taper.*(0.5+0.5*confidence).*double(state.valid);
    [x,y] = meshgrid((state.cols-floor(size(field,2)/2)-1)*inputData.pixelSize, ...
        (state.rows-floor(size(field,1)/2)-1)*inputData.pixelSize);
    radius = sqrt((x-g.M*inputData.pixelSize).^2 + ...
        (y-g.N*inputData.pixelSize).^2 + (g.Z-z)^2);
    state.illumination = cast(white*g.Z/(g.Z-z)* ...
        exp(1i*2*pi/inputData.wavelength*radius),'like',field);
    state.distance = z;
    state.whiteAmplitude = white;
end
u = propagateAngularSpectrum(field,inputData.pixelSize,inputData.wavelength,-state.distance);
t = u(state.rows,state.cols)./state.illumination;
a = abs(t); p = angle(t);
% Morphology is evaluated on CPU; projection retains the field's numeric type.
masks.absorption = morphologySupport(double(gatherIfNeeded(1-a)),state.valid);
masks.phase = morphologySupport(double(gatherIfNeeded(p)),state.valid);
masks.updateWeight = state.updateWeight;
masks.valid = state.valid;
s1 = cast(masks.absorption,'like',a); s2 = cast(masks.phase,'like',p);
alpha = settings.strength.*cast(state.updateWeight,'like',a);
amplitude = a + alpha.*((1-(1-a).*s1)-a);
phase = p.*(1-alpha.*(1-s2));
updated = amplitude.*exp(1i*phase).*state.illumination;
% Preserve unsupported corners and zero-weight boundary pixels exactly in
% the sample plane, avoiding divide/multiply roundoff outside the update.
patch = u(state.rows,state.cols);
active = state.updateWeight>0;
patch(active) = updated(active);
u(state.rows,state.cols) = patch;
field = propagateAngularSpectrum(u,inputData.pixelSize,inputData.wavelength,state.distance);
end

function mask = morphologySupport(values,valid)
% Contrast adjustment, Poisson mixture threshold, alternating morphology,
% then Gaussian smoothing and the automatically selected Sobel threshold.
if any(~isfinite(values(valid)))
    error('Lensless:InvalidConstraintField','Constraint input must be finite.');
end
if max(values(valid)) == min(values(valid))
    % No spatial evidence for a support: preserve this component unchanged.
    mask = true(size(values)); return;
end
% Unsupported corners must not influence thresholds. Nearest valid values
% extend morphology across the exterior without creating a dark aperture.
if any(~valid(:))
    [~,nearest] = bwdist(valid);
    values(~valid) = values(nearest(~valid));
end
values = imadjust(values,stretchlim(values(valid)));
mask = values > poissonThreshold(values(valid));
disk = strel('disk',1); square = strel('square',2);
% Explicit operations retain the finite-array edge convention of the
% experimental implementation (imclose adds different boundary padding).
closed = imerode(imdilate(mask,disk),disk);
opened = imdilate(imerode(mask,disk),disk);
filtered = (double(imdilate(imerode(closed,square),square)) + ...
    double(imerode(imdilate(opened,square),square)))/2;
[~, threshold] = edge(filtered,'sobel');
mask = imbinarize(imgaussfilt(filtered,0.5),threshold);
mask(~valid) = true; % Preserve unsupported components; update weight is zero.
end

function threshold = poissonThreshold(values)
% Minimize the two-class Poisson mixture cost using cumulative histogram
% statistics. Empty classes are excluded; 0*log(0) has limiting value zero.
edges = linspace(0,1,128);
probability = histcounts(mat2gray(values),edges,'Normalization','probability');
weight = cumsum(probability);
moment = cumsum((0:numel(probability)-1).*probability);
other = 1-weight;
valid = weight>0 & other>0;
cost = inf(size(weight));
mu0 = moment(valid)./weight(valid);
mu1 = max(moment(end)-moment(valid),0)./other(valid);
cost(valid) = moment(end)-weight(valid).*log(weight(valid))- ...
    other(valid).*log(other(valid))-weight(valid).*xlogx(mu0)-other(valid).*xlogx(mu1);
[~, index] = min(cost);
threshold = min(values(:))+edges(index)*(max(values(:))-min(values(:)));
end

function result = xlogx(values)
result = zeros(size(values));
positive = values>0;
result(positive) = values(positive).*log(values(positive));
end
