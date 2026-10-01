function verify_simulation_image_bounds()
setup_project;
out=fullfile(tempdir,'lensless_image_bounds');if ~isfolder(out),mkdir(out);end
file=fullfile(out,'rectangle.png');imwrite(uint8(128*ones(20,40)),file);
s=struct('type','image','amplitudeFile',file,'phaseFile','', ...
    'sizePixels',64,'offsetPixels',[0,0],'phaseScale',1,'absorptionScale',1);
o=resolveSimulationCamera(struct('imageSize',100,'sample',s));
assert(o.imageSize==32 && isequal(o.sampleImageSize,[32 64]));
[t,mask]=createSimulationSample(32,48,s);
assert(nnz(mask)==32*64 && all(t(~mask)==0));
assert(all(abs(t(mask)-128/255)<1e-12));
support=false(40);support(9:24,9:24)=true;
m=simulationCameraMasks(support,[8.5 8.5 24.5 24.5;20.5 8.5 36.5 24.5],16);
assert(all(m{1}(:)) && nnz(m{2})==16*4);
% A large requested camera must be clamped BEFORE sampling/grid checks.
o=struct('imageSize',2000,'padding',320,'pixelSize',3e-6,'wavelength',514e-9, ...
    'sourceToSample',30e-3,'sourceOffsetPixels',[300,-180], ...
    'sampleDistances',(.5:1:12.5)*1e-3,'calibrationDistances',(.5:.25:1.5)*1e-3, ...
    'calibration',struct('imageSize',192),'sample',s);
o.sample.amplitudeFile='cameraman.tif';o.sample.sizePixels=512;
o=resolveSimulationCamera(o);assert(o.imageSize==512);
d=checkSimulationSampling(o);assert(all(d.maxCarrierCyclesPerPixel<.5));
o.sourceOffsetPixels=[1000,-28];o.sourceToSample=15e-3;failed=false;
try,checkSimulationSampling(o);catch e,failed=strcmp(e.identifier,'Lensless:SimulationAliasing');end
assert(failed);
fprintf('Camera clamping, finite image support and shifted out-of-image mask checks passed.\n');
end
