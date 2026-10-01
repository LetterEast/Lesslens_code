function options = defaultReconstructionOptions(profile)
%DEFAULTRECONSTRUCTIONOPTIONS Shared APRW settings with explicit profile differences.
% All distances are in metres. Fast retains its historical focus settings.
if nargin < 1, profile = 'standard'; end
profile = validatestring(profile, {'standard','fast'});
projectRoot = fileparts(fileparts(fileparts(mfilename('fullpath'))));
options.iterations = 18;
options.recordEvery = 3;
options.adaptiveConstraint.enabled = false; % Experimental object-plane constraint
options.adaptiveConstraint.strength = 0.05; % 0 = bypass; 1 = full local update (also weighted spatially)
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
    options.showFigures = false;
end
end
