function options = defaultReconstructionOptions(profile)
%DEFAULTRECONSTRUCTIONOPTIONS Shared APRW settings with explicit profile differences.
% All distances are in metres. Fast retains its historical TV/focus settings.
if nargin < 1, profile = 'standard'; end
profile = validatestring(profile, {'standard','fast'});
projectRoot = fileparts(fileparts(fileparts(mfilename('fullpath'))));
options.iterations = 18;
options.recordEvery = 3;
options.adaptiveConstraint.enabled = false; % Experimental object-plane constraint
options.adaptiveConstraint.strength = 0.05; % 0 = bypass; 1 = full local update (also weighted spatially)
options.adaptiveConstraint.phaseMode = 'wrapped'; % Original; 'circular' avoids wrap seams; 'off' amplitude only
options.adaptiveConstraint.edgeWidth = 32; % Sample-grid pixels; taper update at union boundary
options.adaptiveConstraint.distance = []; % [m]; [] autofocus once before first application
options.adaptiveConstraint.whiteAmplitude = []; % [] estimate from first image, 99.9 percentile
options.focus.prior = 1.68e-3;
options.focus.halfRange = 3e-3;
options.focus.step = 0.01e-3;
options.focus.method = 'fft'; % 'fft' (original) or 'adfrft' (fractional domain)
options.focus.frftOrder = 0.8; % ADFrFT order: 0 < p <= 1
options.focus.cropSize = 256; % ADFrFT central ROI side [pixels]; 0 = full field
options.focus.roiMode = 'center'; % 'center', 'full', or interactive 'manual'
options.focus.roiBounds = []; % Optional [xmin ymin xmax ymax] to reuse a manual ROI
options.tv.enabled = true;
options.tv.lambdaMin = 2e-3;  % maximum-overlap protected region
options.tv.lambdaMax = 2e-2;  % single-measurement region
options.tv.coveragePower = 2;
options.tv.coverageWeight = 0.8;
options.tv.uncertaintyWeight = 0.2;
options.tv.boundaryWeight = 0.5;
options.tv.boundaryPower = 2;
options.tv.directionGain = 0.4;
options.tv.step = 2;
options.tv.subiterations = 10;
options.output.rootFolder = fullfile(projectRoot, 'outputs', 'reconstruction');
options.output.cropToValidFOV = true;
options.output.zeroFillInvalid = true;
% Output uses the geometric sample-plane union; unsupported corners are zero-filled.
% originalFOV_amplitude.png uses the first camera footprint mapped to the sample.
options.output.defaultToOriginalFOV = false;
% Retain every geometrically supported sample pixel after back propagation.
options.output.minimumCoverageCount = 1;
options.output.limitHorizontalToReferenceFOV = false;
options.output.trimPropagationBoundary = false;
options.showFigures = true;

options.output.saveFocusPlot = true;
if strcmp(profile, 'fast')
    options.recordEvery = options.iterations;
    options.focus.halfRange = 0.6e-3;
    options.tv.coveragePower = 1;
    options.tv.coverageWeight = 1;
    options.tv.subiterations = 50;
    options.showFigures = false;
end
end
