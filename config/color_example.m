% Color reconstruction reads a manifest of prepared per-channel MAT files.
config.inputFile = fullfile(projectRoot, 'data', 'color_input.mat');
config.options = defaultReconstructionOptions('standard');
% Dataset settings. All lengths are in metres unless stated otherwise.
% Public color configuration: edit this file, then run reconstruct_color[_fast].
% Loaded by loadProjectConfig('color',false); do not run this script directly.
% Personal paths belong in root color_config.m, which is excluded from Git.
%% Dataset configuration
config.outputRoot = fullfile(projectRoot, 'outputs', 'color');
% Set COLOR_FORCE_RECONSTRUCT=1 before launching MATLAB to ignore the cache.
config.reuseExistingReconstructions = ...
    ~strcmp(getenv('COLOR_FORCE_RECONSTRUCT'), '1');
% The input manifest lists channels in B, G, R order.
config.channelNames = {'B', 'G', 'R'};
config.rgbOrder = [3, 2, 1];

%% Reconstruction: independent color settings, shared grayscale/core algorithm
config.options.iterations = 18;
config.options.recordEvery = 3;
config.options.showFigures = true;
config.options.output.saveFocusPlot = true;

%% Focus (metres); each channel searches independently
config.options.focus.method = 'fft'; % 'fft' / 'adfrft'
config.options.focus.prior = 1.68e-3;
config.options.focus.halfRange = 0.6e-3; % 0: fixed distance
config.options.focus.step = 0.01e-3;
config.options.focus.frftOrder = 0.8;
config.options.focus.cropSize = 256;
config.options.focus.roiMode = 'center'; % 'manual' / 'center' / 'full'
config.options.focus.roiBounds = [];

%% TV and adaptive support
config.options.tv.enabled = true;
config.options.tv.lambdaMin = 2e-3;
config.options.tv.lambdaMax = 2e-2;
config.options.adaptiveConstraint.enabled = false;
config.options.adaptiveConstraint.strength = 0.05;
config.options.adaptiveConstraint.phaseMode = 'wrapped'; % 'circular' / 'off' also supported
config.options.adaptiveConstraint.edgeWidth = 32;
config.options.adaptiveConstraint.distance = [];
config.options.adaptiveConstraint.whiteAmplitude = [];

%% Grayscale export; full channel canvas is also kept for fusion
config.options.output.cropToValidFOV = true;
config.options.output.zeroFillInvalid = true;
config.options.output.defaultToOriginalFOV = false;

%% Color fusion

% Correct only broad channel-dependent shading during RGB fusion. The
% reconstructed grayscale channels themselves are left unchanged.
config.colorBalance.enabled = true;
config.colorBalance.backgroundSigma = 100;
config.colorBalance.maximumGain = 1.25;
config.colorBalance.flattenLowCoverage = true;
config.colorBalance.fullColorCoverageFraction = 0.60;

% Second primary output: convert reconstructed field amplitude to
% transmitted intensity and apply one global white gain per wavelength.
% No coverage-dependent colour taper or local spatial colour correction is
% applied to this version.
config.physicalIntensity.whiteReferencePercentile = 85;
config.physicalIntensity.minimumWhitePixels = 2048;
config.physicalIntensity.minimumCoverageFraction = 0.60;
config.physicalIntensity.targetWhite = 0.98;

% Correct scalar exposure drift between measurements before APRW. This is
% not image averaging: every frame is retained, and only one gain is
% applied to each frame. The gain limit prevents a damaged frame from
% dominating the phase-retrieval constraint.
config.inputBrightness.enabled = false; % Match grayscale measurements; enable only for known exposure drift
config.inputBrightness.percentiles = [5, 95];
config.inputBrightness.maximumGain = 3;

%% Same reconstruction and grayscale display defaults as single-channel runs.
% reconstructColor retains full fields separately for RGB registration.
config.options.inputBrightness = config.inputBrightness;

