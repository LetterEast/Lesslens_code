function [result,folder]=simulateModelComparison(options)
%SIMULATEMODELCOMPARISON End-to-end diffraction/calibration/model ablation.
% Calibration and specimen are DIFFERENT synthetic objects. Calibrated
% source parameters come only from grid intensities, not ground-truth centres.
% All arms see identical specimen images. The first three use simplified
% amplitude projection; the optional fourth invokes the experimental core.
% Complete options are supplied by main/demo_sim.m.
o=resolveSimulationCamera(options);rng(o.seed);n=o.imageSize;p=o.padding;side=n+2*p;
sampling=checkSimulationSampling(o);
validateattributes(o.iterations,{'numeric'},{'scalar','integer','positive'});
validateattributes(o.recordEvery,{'numeric'},{'scalar','integer','positive'});
validateattributes(o.sample.phaseScale,{'numeric'},{'scalar','finite','real','nonnegative'});
validateattributes(o.sample.absorptionScale,{'numeric'},{'scalar','finite','real','nonnegative'});
validateattributes(o.evaluation.phaseAmplitudeThreshold,{'numeric'},{'scalar','finite','real','nonnegative'});
fprintf('1/3: Generate independent calibration-grid diffraction images and estimate M/N/Z.\n');
crop=p+(1:n);[x,y]=meshgrid(((1:side)-floor(side/2)-1)*o.pixelSize);
xx=x/o.pixelSize;yy=y/o.pixelSize;
distances=o.sampleDistances(:).';positions=distances-distances(1);
assert(numel(distances)>=3 && all(diff(distances)>0));
source=struct('M',o.sourceOffsetPixels(1),'N',o.sourceOffsetPixels(2), ...
    'Z',o.sourceToSample+distances(1));
H=cell(size(distances));
freq=ifftshift((-floor(side/2):ceil(side/2)-1)/(side*o.pixelSize));
[fx,fy]=meshgrid(freq);root=sqrt(max(1-(o.wavelength*fx).^2-(o.wavelength*fy).^2,0));
for k=1:numel(H),H{k}=exp(1i*2*pi/o.wavelength*distances(k)*root);end
trueModel=makeModel(source,distances(1),o,x,y,distances,H,false);
% Compact absorbing dot grid, with a stronger centre to retain grid labels.
pitch=o.calibration.gridPitchPixels;gridObject=ones(side);
for ix=-2:2
    for iy=-2:2
        radius=4.5;depth=0.65;
        if ix==0 && iy==0,radius=7;depth=0.9;end
        dot=exp(-((xx-ix*pitch).^2+(yy-iy*pitch).^2).^2/(2*radius^4));
        gridObject=gridObject.*(1-depth*dot);
    end
end
calDistances=o.calibrationDistances(:).';
assert(calDistances(1)==distances(1),'Calibration and specimen share the first detector position.');
calH=cell(size(calDistances));
for k=1:numel(calH),calH{k}=exp(1i*2*pi/o.wavelength*calDistances(k)*root);end
calModel=makeModel(source,distances(1),o,x,y,calDistances,calH,false);
calibrationImages=cell(size(calH));
calSize=o.calibration.imageSize;
validateattributes(calSize,{'numeric'},{'scalar','integer','positive','<=',side});
calCrop=floor(side/2)-floor(calSize/2)+(1:calSize);
for k=1:numel(calH)
    u=forward(gridObject,calModel,k);calibrationImages{k}=abs(u(calCrop,calCrop)).^2;
end
calOptions=defaultCalibrationOptions();calOptions.pitchPixels=o.calibration.approxPitchPixels;
calOptions.backgroundSigma=o.calibration.backgroundSigma;
calOptions.projectionSmooth=o.calibration.projectionSmooth;calOptions.showFigures=o.showFigures;
[estimated,calibration]=calibratePointSource(calibrationImages,calDistances-calDistances(1),calOptions);
% Simulation alone knows the ground truth: reject gross calibration failure
% instead of drawing an apparently plausible FOV using a wrong estimate.
if abs(estimated.Z/source.Z-1)>.2 || norm([estimated.M-source.M,estimated.N-source.N])>.2*calSize
    error('Lensless:SimulationCalibrationFailed', ...
        ['Calibration failed: true M/N/Z = %.2f / %.2f / %.3f mm; estimated = %.2f / %.2f / %.3f mm. ' ...
        'Inspect calibration grid visibility and pitch.'],source.M,source.N,source.Z*1e3,estimated.M,estimated.N,estimated.Z*1e3);
end
% Reuse experimental sample-plane geometry. These raw detector frames have
% not been registered, so orig_M/orig_N are zero (no registration shear).
g=struct('M',source.M,'N',source.N,'Z',source.Z,'orig_M',0,'orig_N',0, ...
    'padSize',[p p],'orig_size',[n n]);
trueSupport=sampleGeometrySupport([side side],g,positions,distances(1));
g.M=estimated.M;g.N=estimated.N;g.Z=estimated.Z;
estimatedSupport=sampleGeometrySupport([side side],g,positions,distances(1));
fov=struct('truth',trueSupport,'estimated',estimatedSupport);
fov.addedMask=trueSupport.mask & ~trueSupport.firstMask;
fov.areaGain=nnz(trueSupport.mask)/nnz(trueSupport.firstMask);
fov.estimatedAreaGain=nnz(estimatedSupport.mask)/nnz(estimatedSupport.firstMask);
fov.pixelAreaUm2=(o.pixelSize*1e6)^2;
fprintf('Sample FOV: first %d pixels, union %d pixels, added %.2f%% (estimated %.2f%%).\n', ...
    nnz(trueSupport.firstMask),nnz(trueSupport.mask),100*(fov.areaGain-1),100*(fov.estimatedAreaGain-1));
displayMask=estimatedSupport.firstMask;
if o.fov.displayUnion,displayMask=estimatedSupport.mask;end
[dr,dc]=find(displayMask);displayRows=min(dr):max(dr);displayCols=min(dc):max(dc);
fprintf('2/3: Generate specimen images using the TRUE off-axis spherical illumination.\n');
[truth,sampleSupport]=createSimulationSample(n,p,o.sample);
cameraMasks=simulationCameraMasks(sampleSupport,trueSupport.planeEdges,n);
finiteImage=strcmp(o.sample.type,'image');
if finiteImage
    % A finite image is a zero-extended object: propagate t*illumination
    % directly. Calibration keeps its separate transparent-background model.
    trueModel.referenceValue=0;
    for k=1:numel(H),trueModel.background{k}=zeros(side);end
end
images=cell(size(H));
for k=1:numel(H)
    u=forward(truth,trueModel,k);im=abs(u(crop,crop)).^2;
    images{k}=max(im+o.noiseStd*randn(size(im)),0);
    images{k}(~cameraMasks{k})=0;
end
% Plane baseline gets the SAME estimated lateral source direction. It uses
% a plane wave plus geometric translation, with no curvature/magnification.
models={trueModel,makeModel(estimated,distances(1),o,x,y,distances,H,false), ...
    makeModel(estimated,distances(1),o,x,y,distances,H,true)};
if finiteImage
    for j=2:3
        models{j}.referenceValue=0;
        for k=1:numel(H),models{j}.background{k}=zeros(side);end
    end
end
names={'known_spherical','calibrated_spherical','plane_plus_shift'};
labels={'已知真值球面波（参考组）','标定后球面波（使用估计 M/N/Z）','平面波 + 平移（对照组）'};
if o.core.enabled
    names{end+1}='experimental_core';labels{end+1}='实验核心：APRW + 所选约束';
    models{end+1}=models{2}; % Common estimated model for measurement residual scoring.
end
modelCount=numel(names);coreInfo=struct();
folder=fullfile(o.outputRoot,['model_comparison_' char(datetime('now','Format','yyyyMMdd_HHmmss_SSS'))]);mkdir(folder);
fprintf('3/3: Reconstruct identical specimen measurements with %d comparison arms.\n',modelCount);
fields=cell(1,modelCount);rawFields=cell(1,modelCount);phaseOffsets=zeros(1,modelCount);
histories=zeros(o.iterations,modelCount);metrics=zeros(modelCount,3);
evaluation=trueSupport.mask & sampleSupport;
regionNames={'first_fov','union_fov','added_fov','overlap_fov'};
regionMasks={trueSupport.firstMask & sampleSupport,evaluation, ...
    fov.addedMask & sampleSupport,(trueSupport.coverage==numel(H)) & sampleSupport};
fov.sampleUnionMask=evaluation;fov.sampleFirstMask=regionMasks{1};
fov.sampleAreaGain=nnz(evaluation)/max(nnz(regionMasks{1}),1);
qualityRows=struct([]);
monitor=[];if o.showFigures,monitor=figure('Name','Simulation: model comparison');end
timer=tic;
for modelIndex=1:modelCount
    m=models{modelIndex};t=double(sampleSupport);
    if strcmp(names{modelIndex},'experimental_core')
        [t,coreInfo]=runCoreSimulationComparison(images,cameraMasks,estimated,o,folder);
        histories(:,modelIndex)=coreInfo.rHistory;
    else
    for iteration=1:o.iterations
        proposal=zeros(side);errorSum=0;measuredSum=0;
        for k=1:numel(H)
            u=forward(t,m,k);a=sqrt(images{k});region=u(crop,crop);
            valid=cameraMasks{k};
            errorSum=errorSum+sum((abs(region(valid))-a(valid)).^2);measuredSum=measuredSum+sum(a(valid).^2);
            region(valid)=a(valid).*exp(1i*angle(region(valid)));u(crop,crop)=region;
            if m.plane
                correction=ifft2(fft2(u-m.background{k})./m.transfer{k});
            else
                correction=ifft2(fft2(u-m.background{k}).*conj(H{k}));
            end
            proposal=proposal+m.referenceValue+correction./m.illumination;
        end
        assert(measuredSum>0,'No nonzero measurements inside the sample. Check image and camera positions.');
        t=proposal/numel(H);t(~sampleSupport)=0;
        histories(iteration,modelIndex)=sqrt(errorSum/measuredSum);
        if o.showFigures && isgraphics(monitor) && (mod(iteration,o.recordEvery)==0 || iteration==1)
            figure(monitor);subplot(1,3,1);showMasked(abs(t),displayMask,displayRows,displayCols,[0 1.2]);colormap(gca,gray(256));title(labels{modelIndex},'FontName','Microsoft YaHei');
            subplot(1,3,2);showMasked(angle(t),displayMask,displayRows,displayCols,[-1 1]);title(sprintf('Iteration %d',iteration));
            subplot(1,3,3);plot(histories(1:iteration,modelIndex));grid on;title('Amplitude residual');drawnow;
        end
    end
    end
    % One global phase alignment per model, shared by ALL regions. This is
    % evaluation-only and uses truth, never fed back into reconstruction.
    rawFields{modelIndex}=t;
    phaseValid=evaluation & abs(truth)>o.evaluation.phaseAmplitudeThreshold;
    offset=angle(sum(t(phaseValid).*conj(truth(phaseValid))));t=t*exp(-1i*offset);
    t(~sampleSupport)=0; % Output rule only for the unchanged experimental core.
    phaseOffsets(modelIndex)=offset;
    fields{modelIndex}=t;
    metrics(modelIndex,1)=sqrt(mean((abs(t(evaluation))-abs(truth(evaluation))).^2));
    metrics(modelIndex,2)=sqrt(mean(angle(t(phaseValid).*conj(truth(phaseValid))).^2));
    % Final-field residual (histories contain pre-update residuals).
    errorSum=0;measuredSum=0;
    for k=1:numel(H)
        u=forward(rawFields{modelIndex},m,k);a=sqrt(images{k});
        region=u(crop,crop);valid=cameraMasks{k};
        errorSum=errorSum+sum((abs(region(valid))-a(valid)).^2);measuredSum=measuredSum+sum(a(valid).^2);
    end
    metrics(modelIndex,3)=sqrt(errorSum/measuredSum);
    for j=1:numel(regionNames)
        scores=evaluateSimulationField(t,truth,regionMasks{j},o.evaluation.phaseAmplitudeThreshold);
        row=struct('model',names{modelIndex},'region',regionNames{j});
        keys=fieldnames(scores);for q=1:numel(keys),row.(keys{q})=scores.(keys{q});end
        row.amplitude_residual=metrics(modelIndex,3);
        if isempty(qualityRows),qualityRows=row;
        else,qualityRows(end+1)=row;end %#ok<AGROW>
    end
    fprintf('%s: amplitude RMSE %.5g, wrapped phase RMSE %.5g rad\n',names{modelIndex},metrics(modelIndex,1:2));
end
elapsed=toc(timer);
quality=struct2table(qualityRows);
disp(quality(:,{'model','region','amplitude_psnr_db','amplitude_ssim','phase_rmse_rad'}));
result=struct('truth',truth,'fields',{fields},'names',{names},'labels',{labels},'metrics',metrics, ...
    'coreInfo',coreInfo, ...
    'rawFields',{rawFields},'phaseOffsets',phaseOffsets, ...
    'histories',histories,'estimatedSource',estimated,'trueSource',source, ...
    'calibration',calibration,'images',{images},'calibrationImages',{calibrationImages}, ...
    'crop',crop,'evaluationMask',evaluation,'fov',fov,'quality',quality,'sampling',sampling, ...
    'sampleSupport',sampleSupport,'cameraMasks',{cameraMasks}, ...
    'displayRows',displayRows,'displayCols',displayCols,'options',o,'elapsedSeconds',elapsed);
save(fullfile(folder,'simulation.mat'),'-struct','result');
tableOut=array2table(metrics,'VariableNames',{'amplitude_rmse','phase_rmse_rad','amplitude_residual'},'RowNames',names);
unionScores=quality(strcmp(quality.region,'union_fov'),:);
extra={'amplitude_psnr_db','amplitude_ssim','phase_psnr_db','phase_phasor_ssim','complex_nrmse'};
for k=1:numel(extra),tableOut.(extra{k})=unionScores.(extra{k});end
writetable(tableOut,fullfile(folder,'metrics.csv'),'WriteRowNames',true);
writetable(quality,fullfile(folder,'quality_by_region.csv'));
fovTable=table({'truth';'calibrated'}, ...
    [nnz(trueSupport.firstMask);nnz(estimatedSupport.firstMask)], ...
    [nnz(trueSupport.mask);nnz(estimatedSupport.mask)], ...
    [fov.areaGain;fov.estimatedAreaGain], ...
    'VariableNames',{'geometry','first_fov_pixels','union_fov_pixels','area_gain'});
fovTable.sample_first_pixels=[nnz(trueSupport.firstMask & sampleSupport);nnz(estimatedSupport.firstMask & sampleSupport)];
fovTable.sample_union_pixels=[nnz(trueSupport.mask & sampleSupport);nnz(estimatedSupport.mask & sampleSupport)];
fovTable.sample_area_gain=fovTable.sample_union_pixels./max(fovTable.sample_first_pixels,1);
writetable(fovTable,fullfile(folder,'fov_summary.csv'));
parameters={'M_pixels';'N_pixels';'Z_mm'};
trueValues=[source.M;source.N;source.Z*1e3];estimatedValues=[estimated.M;estimated.N;estimated.Z*1e3];
writetable(table(parameters,trueValues,estimatedValues,estimatedValues-trueValues, ...
    'VariableNames',{'parameter','truth','calibrated','error'}),fullfile(folder,'calibration_parameters.csv'));
% The main comparison remains visible after a visual run. In particular,
% do not leave only the last (plane-wave) iteration monitor on screen.
visible='off';if o.showFigures,visible='on';end
f=figure('Visible',visible,'Name','仿真结果：标定球面波与平面波对比','Position',[80 80 1350 680]);
displayFields={truth,fields{2},fields{3}};
displayLabels={'真实样品（幅相真值）',labels{2},labels{3}};
for k=1:3
    a=displayFields{k};subplot(2,3,k);showMasked(abs(a),displayMask,displayRows,displayCols,[0 1.2]);
    colormap(gca,gray(256));
    label=displayLabels{k};
    if k>1
        q=quality(strcmp(quality.model,names{k}) & strcmp(quality.region,'union_fov'),:);
        label=sprintf('%s\nUnion PSNR %.2f dB | SSIM %.3f',label,q.amplitude_psnr_db,q.amplitude_ssim);
    end
    title(label,'FontName','Microsoft YaHei');
    subplot(2,3,k+3);showMasked(angle(a),displayMask,displayRows,displayCols,[-1 1]);colorbar;
    title('相位 [rad]','FontName','Microsoft YaHei');
end
sgtitle(sprintf('同一组离焦图，%d 次迭代；标定估计 M=%.2f px，N=%.2f px，Z=%.3f mm', ...
    o.iterations,estimated.M,estimated.N,estimated.Z*1e3),'FontName','Microsoft YaHei');
exportgraphics(f,fullfile(folder,'01_calibrated_vs_plane.png'),'Resolution',130);
if ~o.showFigures,close(f);end
f=figure('Visible','off','Position',[100 100 1300 620]);
allFields=[{truth},fields];titles=[{'真实样品'},labels];columns=numel(allFields);
for k=1:columns
    a=allFields{k};subplot(2,columns,k);showMasked(abs(a),displayMask,displayRows,displayCols,[0 1.2]);colormap(gca,gray(256));title(titles{k},'Interpreter','none','FontName','Microsoft YaHei');
    subplot(2,columns,k+columns);showMasked(angle(a),displayMask,displayRows,displayCols,[-1 1]);title('Phase [rad]');colorbar;
end
sgtitle('Amplitude [0, 1.2]; phase [-1, 1] rad; identical data and solver');
exportgraphics(f,fullfile(folder,'amplitude_phase_comparison.png'),'Resolution',130);close(f);
plotSimulationFOV(result,folder);
plotSimulationAcquisition(result,folder);
if o.core.enabled,plotCoreSimulationComparison(result,folder);end
f=figure('Visible','off','Position',[100 100 1000 540]);
for k=1:numel(images)
    subplot(2,ceil(numel(images)/2),k);imshow(images{k},[0 1.2]);
    title(sprintf('z = %.2f mm',distances(k)*1e3));
end
exportgraphics(f,fullfile(folder,'simulated_measurements.png'),'Resolution',130);close(f);
f=figure('Visible','off');
subplot(1,2,1);plot(calibration.positions*1e3,calibration.pitchPixels,'o');hold on;
plot(calibration.positions*1e3,polyval(calibration.pitchFit,calibration.positions));
xlabel('Detector offset [mm]');ylabel('Grid pitch [pixels]');grid on;
subplot(1,2,2);imshow(calibrationImages{1},[]);title('Independent calibration grid');
exportgraphics(f,fullfile(folder,'calibration_fit.png'),'Resolution',130);close(f);
f=figure('Visible','off');hold on;
for k=1:modelCount,plot(1:o.iterations,histories(:,k),'LineWidth',1.3);end
legend(labels,'Interpreter','none','FontName','Microsoft YaHei');grid on;
xlabel('Iteration');ylabel('Internal amplitude residual');
title('Solver-specific histories; compare final residuals in metrics.csv');
exportgraphics(f,fullfile(folder,'convergence.png'));close(f);
fprintf('Calibration truth/estimate: M %.3f / %.3f px, N %.3f / %.3f px, Z %.4f / %.4f mm\n', ...
    source.M,estimated.M,source.N,estimated.N,source.Z*1e3,estimated.Z*1e3);
fprintf('Simulation results: %s\n',folder);
fprintf('Open first: %s\n',fullfile(folder,'01_calibrated_vs_plane.png'));
if o.core.enabled,fprintf('Core comparison: %s\n',fullfile(folder,'06_core_vs_calibrated.png'));end
end

function showMasked(values,mask,rows,cols,limits)
im=imagesc(values(rows,cols),limits);im.AlphaData=double(mask(rows,cols));
axis image;set(gca,'Color',[.7 .7 .7],'XTick',[],'YTick',[],'Box','off');
end

function model=makeModel(source,firstDistance,o,x,y,distances,H,plane)
L=source.Z-firstDistance;assert(L>0,'Estimated source must lie behind the specimen.');
sx=source.M*o.pixelSize;sy=source.N*o.pixelSize;k0=2*pi/o.wavelength;
model=struct('plane',plane,'referenceValue',1,'background',{{}},'transfer',{{}});
radius=sqrt((x-sx).^2+(y-sy).^2+L^2);
model.illumination=L./radius.*exp(1i*k0*(radius-L));
side=size(x,1);freq=ifftshift((-floor(side/2):ceil(side/2)-1)/(side*o.pixelSize));[fx,fy]=meshgrid(freq);
if plane,model.illumination=ones(size(x));end
for j=1:numel(distances)
    d=distances(j);
    if plane
        gain=L/(L+d);shiftX=-d*sx/L;shiftY=-d*sy/L;
        model.transfer{j}=gain*H{j}.*exp(-1i*2*pi*(fx*shiftX+fy*shiftY));
        model.background{j}=gain*ones(size(x));
    else
        r=sqrt((x-sx).^2+(y-sy).^2+(L+d)^2);
        model.background{j}=L./r.*exp(1i*k0*(r-L));
        model.transfer{j}=H{j};
    end
end
end

function u=forward(t,model,k)
if model.plane
    u=model.background{k}+ifft2(fft2(t-model.referenceValue).*model.transfer{k});
else
    % transfer is set by caller below through the distance-specific kernel.
    u=model.background{k}+ifft2(fft2((t-model.referenceValue).*model.illumination).*model.transfer{k});
end
end
