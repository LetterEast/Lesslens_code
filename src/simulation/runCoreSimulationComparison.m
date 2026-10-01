function [transmission,info]=runCoreSimulationComparison(images,masks,source,o,folder)
%RUNCORESIMULATIONCOMPARISON Adapt raw simulation frames to the REAL core API.
% Frames already share one laboratory coordinate grid: no registration shear.
% M/N stay nonzero; orig_M/orig_N=0 prevents adding a fictitious translation.
n=o.imageSize;p=o.padding;side=n+2*p;crop=p+(1:n);
g=struct('M',source.M,'N',source.N,'Z',source.Z,'orig_M',0,'orig_N',0, ...
    'orig_size',[n n],'padSize',[p p]);
inputData.images=cell(size(images));g.ValidMaskHard=cell(size(images));
for k=1:numel(images)
    inputData.images{k}=zeros(side);inputData.images{k}(crop,crop)=images{k};
    g.ValidMaskHard{k}=false(side);g.ValidMaskHard{k}(crop,crop)=masks{k};
end
inputData.geometry=g;inputData.distanceSteps=[0,diff(o.sampleDistances)];
inputData.wavelength=o.wavelength;inputData.pixelSize=o.pixelSize;
validateReconstructionInput(inputData);

settings=defaultReconstructionOptions();
settings.iterations=o.iterations;settings.recordEvery=o.iterations;
settings.stopOnConvergence=false; % Same iteration budget as the baseline.
settings.showFigures=o.showFigures;
settings.focus.prior=o.sampleDistances(1);settings.focus.halfRange=0;
settings.focus.roiMode='full';
settings.tv.enabled=o.core.tvEnabled;
settings.adaptiveConstraint.enabled=o.core.adaptiveConstraintEnabled;
settings.adaptiveConstraint.strength=o.core.adaptiveConstraintStrength;
settings.adaptiveConstraint.edgeWidth=o.core.adaptiveConstraintEdgeWidth;
settings.adaptiveConstraint.distance=o.sampleDistances(1);
L=source.Z-o.sampleDistances(1);
% Simulation intensity units: incident amplitude on-axis at plane 1 is L/Z.
settings.adaptiveConstraint.whiteAmplitude=L/source.Z;
settings.output.rootFolder=fullfile(folder,'experimental_core');
settings.output.cropToValidFOV=false;settings.output.zeroFillInvalid=false;
settings.output.saveFocusPlot=false;
timer=tic;[~,runFolder]=reconstructMultiPlane(inputData,settings);elapsed=toc(timer);
manifest=jsondecode(fileread(fullfile(runFolder,'run_manifest.json')));
checkpoint=fullfile(runFolder,manifest.latestResult);
saved=load(checkpoint,'object','focusDistance');
meta=load(fullfile(fileparts(fileparts(checkpoint)),'diagnostics','meta.mat'),'rHistory');
assert(abs(saved.focusDistance-o.sampleDistances(1))<1e-12);
assert(isequal(size(saved.object),[side side]));
% Core saves the sample field after phase-only illumination removal and TV.
% Convert to the SAME dimensionless transmission used by the other arms.
% This is a geometric 1/r correction, not a fitted gain or display scaling.
[x,y]=meshgrid(((1:side)-floor(side/2)-1)*o.pixelSize);
radius=sqrt((x-source.M*o.pixelSize).^2+(y-source.N*o.pixelSize).^2+L^2);
transmission=saved.object./(L./radius);
info=struct('algorithm','reconstructMultiPlane','runFolder',runFolder, ...
    'settings',settings,'iterationsCompleted',manifest.iterationsCompleted, ...
    'rHistory',meta.rHistory,'elapsedSeconds',elapsed, ...
    'coordinateConvention','unregistered laboratory grid, orig_M=orig_N=0', ...
    'normalization','sample output divided by estimated spherical amplitude L/r', ...
    'focusMode','fixed known simulation distance', ...
    'supportDifference','baseline enforces finite image support each iteration; core only output is masked');
save(fullfile(runFolder,'simulation_adapter.mat'),'g','settings');
end
