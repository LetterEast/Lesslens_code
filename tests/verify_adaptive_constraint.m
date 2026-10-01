function verify_adaptive_constraint()
% Public synthetic regression: bypass, run isolation, and APRW integration.
root = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root,'src'));
rng(7);
g = struct('Z',0.05,'M',0,'N',0,'orig_M',0,'orig_N',0, ...
    'padSize',[0,0],'orig_size',[32,40]);
g.ValidMaskHard = {true(32,40)};
d = struct('geometry',g,'images',{{0.2+0.6*rand(32,40)}}, ...
    'distanceSteps',0,'pixelSize',3e-6,'wavelength',514e-9);
o = defaultReconstructionOptions();
s = o.adaptiveConstraint; s.distance = 1e-3; s.strength = 0;
u = rand(32,40).*exp(1i*randn(32,40));
assert(isequal(u,applyAdaptiveConstraint(u,d,s)));
s.strength = 0.05;
[v,state] = applyAdaptiveConstraint(u,d,s);
assert(all(isfinite(v(:))) && norm(v-u,'fro')>0);
assert(isequal(v,applyAdaptiveConstraint(u,d,s,state)));
s.whiteAmplitude = 0.5;
[~,other] = applyAdaptiveConstraint(u,d,s,[]);
assert(other.whiteAmplitude==0.5 && state.whiteAmplitude~=0.5);
s.strength = 1;
assert(all(isfinite(applyAdaptiveConstraint(u,d,s)), 'all'));
o.iterations = 1; o.recordEvery = 1; o.showFigures = false;
o.tv.enabled = false; o.output.saveFocusPlot = false;
o.focus.prior = 1e-3; o.focus.halfRange = 0;
o.output.rootFolder = fullfile(root,'outputs','constraint_verification');
[baseline,~] = reconstructMultiPlane(d,o);
o.adaptiveConstraint.enabled = true; o.adaptiveConstraint.strength = 0;
[bypass,~] = reconstructMultiPlane(d,o);
assert(isequal(baseline,bypass));
o.adaptiveConstraint.strength = 0.05;
[fields,folder] = reconstructMultiPlane(d,o);
assert(all(isfinite(fields{1}(:))) && ~isequal(fields,baseline));
saved = load(fullfile(folder,'constraint_0001.mat'));
focus = load(fullfile(folder,'constraint_autofocus.mat'));
assert(saved.constraintSettings.distance==focus.bestDistance);
assert(saved.constraintSettings.whiteAmplitude>0);
assert(all(isfield(saved.constraintMasks,{'absorption','phase'})));
verifyUnionUpdate(d,o.adaptiveConstraint);
fprintf('Adaptive constraint checks passed.\n');
end

function verifyUnionUpdate(d,s)
% Shifted footprints: old first-camera boundary lies inside the valid union.
g = struct('Z',0.05,'M',0,'N',0,'orig_M',0,'orig_N',-100, ...
    'padSize',[28,20],'orig_size',[40,60]);
d.geometry = g; d.distanceSteps = [0,0.02];
d.images = {ones(96,100),ones(96,100)};
s.distance = 0; s.whiteAmplitude = 1; s.edgeWidth = 8; s.strength = 0.15;
u = rand(96,100).*exp(1i*randn(96,100));
[v,state,masks] = applyAdaptiveConstraint(u,d,s);
support = sampleGeometrySupport(size(u),g,cumsum(d.distanceSteps),0);
weights = zeros(size(u)); weights(state.rows,state.cols) = masks.updateWeight;
assert(any(weights(support.mask & ~support.firstMask)>0));
assert(all(weights(~support.mask)==0));
assert(isequal(u(weights==0),v(weights==0))); % z=0 isolates sample update
assert(all(weights(:)>=0 & weights(:)<=1));
% The first footprint begins at row 29: no on/off seam at this internal edge.
assert(weights(28,50)>0.5 && weights(29,50)>0.5);
assert(abs(weights(29,50)-weights(28,50))<0.1);
assert(all(isfinite(v(:))));
end
