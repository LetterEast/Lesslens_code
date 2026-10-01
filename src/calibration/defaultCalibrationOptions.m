function options=defaultCalibrationOptions()
%DEFAULTCALIBRATIONOPTIONS MainV5 projection/centroid/grid regression settings.
options.pitchPixels=100; % Approximate observed grid spacing at first detector
options.backgroundSigma=100;
options.projectionSmooth=20;
options.showFigures=false;
end
