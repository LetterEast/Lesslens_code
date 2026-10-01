function config=calibration_example()
% Public calibration template; pass its return value to demo_calibration.
root=fileparts(fileparts(mfilename('fullpath')));addpath(root);setup_project;
config.imageFolder=''; % Folder of dark dot-grid images, sorted by filename
config.distanceSteps=[0,1,1,1]*1e-3; % Adjacent detector steps [m], one per image
config.options=defaultCalibrationOptions();
config.options.pitchPixels=100;
config.options.showFigures=true;
config.outputRoot=fullfile(root,'outputs','calibration');
end
