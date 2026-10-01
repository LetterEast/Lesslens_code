function verify_geometry_output()
%VERIFY_GEOMETRY_OUTPUT Regression checks for geometry and APRW output modes.
root = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root,'src'));
g = struct('Z',0.05,'M',0,'N',0,'orig_M',8,'orig_N',-6, ...
    'padSize',[8,10],'orig_size',[16,24]);
imageSize = [32,44]; z = [0,0.001];
identity = sampleGeometrySupport(imageSize,g,0,0);
assert(isequal(identity.bounds,[11,9,34,24]));
assert(nnz(identity.mask)==16*24);
g.ValidMaskHard = {identity.mask,identity.mask};
g.ValidMask = g.ValidMaskHard;
input = struct('images',{{0.25*ones(imageSize),0.25*ones(imageSize)}}, ...
    'geometry',g,'distanceSteps',z,'wavelength',514e-9,'pixelSize',3e-6);
options.iterations = 1; options.recordEvery = 1;
options.focus = struct('prior',0.001,'halfRange',0,'step',1e-5);
options.showFigures = false;
options.output = struct('minimumCoverageCount',1,'trimPropagationBoundary',false);
testRoot = fullfile(root,'outputs','geometry_verification', ...
    char(datetime('now','Format','yyyyMMdd_HHmmss_SSS')));
support = sampleGeometrySupport(imageSize,g,z,0.001);
for mode = 1:3
    options.output.rootFolder = fullfile(testRoot,sprintf('mode%d',mode));
    options.output.cropToValidFOV = mode~=2;
    options.output.zeroFillInvalid = mode~=2;
    options.output.defaultToOriginalFOV = mode==3;
    [~,folder] = reconstructMultiPlane(input,options);
    resultFolder = fullfile(folder,'iteration_0001','results');
    r = load(fullfile(resultFolder,'reconstruction.mat'));
    b = r.bounds;
    expected = support.mask(b(2):b(4),b(1):b(3));
    assert(isequal(r.validMask,expected));
    assert(isequal(r.sampleCoverage,support.coverage(b(2):b(4),b(1):b(3))));
    assert(isequal(size(r.object),size(r.validMask)));
    if mode~=2, assert(isequal(b,support.bounds)); end
    [amplitudePNG,~,alpha] = imread(fullfile(resultFolder,'amplitude.png'));
    [phasePNG,~,phaseAlpha] = imread(fullfile(resultFolder,'phase_heatmap.png'));
    if mode==3, expected = r.originalValidMask; end
    assert(isempty(alpha) && isempty(phaseAlpha));
    assert(all(amplitudePNG(~expected)==0));
    assert(all(phasePNG(repmat(~expected,1,1,3))==0));
    if mode==2, assert(isequal(size(r.object),imageSize)); end
end
fprintf('Geometry output checks passed.\n');
end
