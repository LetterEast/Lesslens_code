function verify_workflow_profiles()
% Verify that visual/fast entry points change display, not physical results.
root=fileparts(fileparts(mfilename('fullpath')));addpath(root);setup_project;
old=get(groot,'defaultFigureVisible');set(groot,'defaultFigureVisible','off');
cleanup=onCleanup(@() set(groot,'defaultFigureVisible',old));
inputData=syntheticInput();out=fullfile(root,'outputs','workflow_verification');
if ~isfolder(out),mkdir(out);end
cfg.inputFile=fullfile(out,'input.mat');save(cfg.inputFile,'inputData');
cfg.options=defaultReconstructionOptions();cfg.options.iterations=1;
cfg.options.recordEvery=1;cfg.options.focus.halfRange=0;cfg.options.tv.enabled=false;
cfg.options.output.rootFolder=out;cfg.options.output.saveFocusPlot=false;
[a,folder]=demo_exp(cfg);[b,~]=demo_exp_fast(cfg);
assert(isequal(a,b),'Visual/fast experiment profiles changed numerical output.');
assert(isfile(fullfile(folder,'run_manifest.json')));
close all;
o=struct('iterations',3,'outputRoot',out,'sourceOffsetPixels',[30,-22]);
o.imageSize=192;o.padding=64;o.sourceToSample=15e-3;
o.sampleDistances=(.5:.5:3)*1e-3;
o.fov=struct('saveAcquisitionGif',false);
o.core=struct('enabled',false);
o.sample=struct('type','mixed'); % Independent of private image selection in the entry.
[a,~]=demo_sim(o);[b,~]=demo_sim_fast(o);
for k=1:3,assert(norm(a.fields{k}-b.fields{k},'fro')<1e-10);end
assert(abs(a.estimatedSource.Z/a.trueSource.Z-1)<0.1);
% Calibration must retain original z positions when a middle frame fails.
images=a.calibrationImages;images{3}=ones(size(images{3}));
[g,diag]=calibratePointSource(images,a.calibration.positions,a.calibration.options);
assert(~diag.validFrames(3));assert(isequal(diag.positions,a.calibration.positions));
assert(isfinite(g.Z) && g.Z>0);
% The folder entry must exercise the same calibration after image import.
imageFolder=fullfile(out,'calibration_images');if ~isfolder(imageFolder),mkdir(imageFolder);end
for k=1:numel(b.calibrationImages)
    image=b.calibrationImages{k};image=image/max(image(:));
    imwrite(uint16(image*65535),fullfile(imageFolder,sprintf('frame_%03d.png',k)));
end
cfgCal=struct('imageFolder',imageFolder,'distanceSteps',[0;diff(b.calibration.positions)], ...
    'options',b.calibration.options,'outputRoot',fullfile(out,'calibration'));
cfgCal.options.showFigures=false;
[g,calFolder]=demo_calibration(cfgCal);
assert(abs(g.Z/b.estimatedSource.Z-1)<0.01);
assert(isfile(fullfile(calFolder,'MNZ_result.mat')));
close all;
fprintf('Visual/fast equivalence and skipped-frame calibration checks passed.\n');
end
