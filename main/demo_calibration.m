function [geometry,folder]=demo_calibration(config)
%DEMO_CALIBRATION Point-source calibration from dark dot-grid images.
% Edit root calibration_config.m once, then run this entry from any directory.
root=fileparts(fileparts(mfilename('fullpath')));addpath(root);setup_project;
if nargin<1
    if ~isfile(fullfile(root,'calibration_config.m'))
        error('Lensless:CalibrationConfig','Set up calibration_config.m, or pass a config based on config/calibration_example.m.');
    end
    config=calibration_config();
end
entries=dir(config.imageFolder);entries=entries(~[entries.isdir]);keep=false(size(entries));
for k=1:numel(entries)
    [~,~,ext]=fileparts(entries(k).name);keep(k)=ismember(lower(ext),{'.png','.jpg','.jpeg','.bmp','.tif','.tiff'});
end
entries=entries(keep);[~,order]=sort(lower(string({entries.name})));entries=entries(order);
if isempty(config.distanceSteps) && isfield(config,'planeSpacing')
    config.distanceSteps=[0,repmat(config.planeSpacing,1,numel(entries)-1)];
end
assert(~isempty(entries),'No calibration images found in the configured folder.');
assert(numel(entries)==numel(config.distanceSteps),'Calibration distanceSteps must match image count.');
positions=cumsum(config.distanceSteps(:));positions=positions-positions(1);
images=cell(1,numel(entries));
for k=1:numel(entries),images{k}=im2double(imread(fullfile(entries(k).folder,entries(k).name)));end
[geometry,diagnostics]=calibratePointSource(images,positions,config.options);
folder=fullfile(config.outputRoot,['calibration_' char(datetime('now','Format','yyyyMMdd_HHmmss_SSS'))]);mkdir(folder);
MNZ_result=geometry; %#ok<NASGU> legacy interchange with prepared measurements
save(fullfile(folder,'calibration.mat'),'geometry','diagnostics','config','entries','MNZ_result');
save(fullfile(folder,'MNZ_result.mat'),'MNZ_result');
f=figure('Visible','off');subplot(1,2,1);plot(positions*1e3,diagnostics.pitchPixels,'o-');
xlabel('Detector offset [mm]');ylabel('Grid pitch [pixels]');grid on;
subplot(1,2,2);plot(diagnostics.centresPixels(:,1),diagnostics.centresPixels(:,2),'o-');
xlabel('Grid centre x [pixels]');ylabel('Grid centre y [pixels]');axis equal;grid on;
exportgraphics(f,fullfile(folder,'calibration_fit.png'));close(f);
fprintf('Source offset M=%.4f px, N=%.4f px; source-to-first-detector Z=%.6f mm\n',geometry.M,geometry.N,geometry.Z*1e3);
fprintf('Calibration saved to: %s\n',folder);
end
