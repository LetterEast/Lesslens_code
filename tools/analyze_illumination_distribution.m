clc;
close all;
clear;
fprintf('Running OL supplement plotting version 2\n');

%% Configuration
dataRoot = "\\192.168.2.166\d\lesslens\2026.9.1\NONE\NONE_G";
imageSubfolder = fullfile("pixel_on", "exposure_-1");
outputFolder = fullfile(fileparts(mfilename('fullpath')), ...
    'ResultFolder', 'BrightnessAnalysis_20260901');
validExtensions = [".png", ".jpg", ".jpeg", ".bmp", ".tif", ".tiff"];
fitDegree = 4;

%% Publication figure style (Optics Letters / Optica supplement)
% Design figures close to their final size in the manuscript so that
% Word/PDF scaling does not make the text too small.
pub.fontName       = 'Arial';
pub.tickFontSize   = 8;      % axis/colorbar tick labels, pt
pub.labelFontSize  = 9;      % x/y/z/colorbar labels, pt
pub.axisLineWidth  = 0.8;    % axes line width, pt
pub.lineWidth      = 1.2;    % curve line width
pub.markerSize     = 4.5;
pub.singleWidthCm  = 9.4;    % compact single-column-like figure width
pub.mediumWidthCm  = 11.5;   % useful when a legend is outside the axes
pub.height2DCm     = 6.8;
pub.height3DCm     = 7.4;
pub.pngResolution  = 600;    % high-resolution raster copy

if ~isfolder(dataRoot)
    error('Data folder does not exist or cannot be accessed: %s', dataRoot);
end
if ~isfolder(outputFolder), mkdir(outputFolder); end

%% Extract brightness from every Pixel_<x>_<y> coordinate
pixelFolders = dir(fullfile(dataRoot, 'Pixel_*_*'));
pixelFolders = pixelFolders([pixelFolders.isdir]);
coordinateX = [];
coordinateY = [];
meanBrightness = [];
maxImageBrightness = [];
imageCount = [];
brightestImage = strings(0,1);
folderName = strings(0,1);

for iFolder = 1:numel(pixelFolders)
    token = regexp(pixelFolders(iFolder).name, ...
        '^Pixel_(-?\d+)_(-?\d+)$', 'tokens', 'once');
    if isempty(token), continue; end

    x = str2double(token{1});
    y = str2double(token{2});
    imageFolder = fullfile(pixelFolders(iFolder).folder, ...
        pixelFolders(iFolder).name, imageSubfolder);
    if ~isfolder(imageFolder)
        warning('Image folder is missing: %s', imageFolder);
        continue;
    end

    files = dir(imageFolder);
    files = files(~[files.isdir]);
    extensions = strings(size(files));
    for iFile = 1:numel(files)
        [~,~,extensions(iFile)] = fileparts(files(iFile).name);
    end
    extensions = lower(extensions);
    files = files(ismember(extensions, validExtensions));
    if isempty(files)
        warning('No supported images found in: %s', imageFolder);
        continue;
    end
    [~, fileOrder] = sort(lower(string({files.name})));
    files = files(fileOrder);

    perImageBrightness = zeros(numel(files), 1);
    for iImage = 1:numel(files)
        imagePath = fullfile(files(iImage).folder, files(iImage).name);
        imageData = im2double(imread(imagePath));
        if size(imageData, 3) == 3, imageData = rgb2gray(imageData); end
        perImageBrightness(iImage) = mean(imageData(:), 'omitnan');
    end

    [thisMaximum, maximumIndex] = max(perImageBrightness);
    coordinateX(end+1,1) = x; %#ok<SAGROW>
    coordinateY(end+1,1) = y; %#ok<SAGROW>
    meanBrightness(end+1,1) = mean(perImageBrightness); %#ok<SAGROW>
    maxImageBrightness(end+1,1) = thisMaximum; %#ok<SAGROW>
    imageCount(end+1,1) = numel(files); %#ok<SAGROW>
    brightestImage(end+1,1) = string(fullfile( ...
        files(maximumIndex).folder, files(maximumIndex).name)); %#ok<SAGROW>
    folderName(end+1,1) = string(pixelFolders(iFolder).name); %#ok<SAGROW>
end

if isempty(meanBrightness)
    error('No images were successfully read under %s.', dataRoot);
end

%% Sort, normalize, find the maximum, and save statistics
[~, pointOrder] = sortrows([coordinateX, coordinateY], [1 2]);
coordinateX = coordinateX(pointOrder);
coordinateY = coordinateY(pointOrder);
meanBrightness = meanBrightness(pointOrder);
maxImageBrightness = maxImageBrightness(pointOrder);
imageCount = imageCount(pointOrder);
brightestImage = brightestImage(pointOrder);
folderName = folderName(pointOrder);

[globalMaximum, globalMaximumIndex] = max(meanBrightness);
if globalMaximum <= 0 || ~isfinite(globalMaximum)
    error('The global mean brightness is invalid: %.8g', globalMaximum);
end
normalizedBrightness = meanBrightness / globalMaximum;

%% Fit a complete 2-D polynomial surface and search integer coordinates
% Center and scale both coordinates to improve numerical conditioning.
% A complete degree-N model contains every u^p*v^q term with p+q <= N.
xCenter = mean(coordinateX);
yCenter = mean(coordinateY);
xScale = max(std(coordinateX), 1);
yScale = max(std(coordinateY), 1);
u = (coordinateX - xCenter) / xScale;
v = (coordinateY - yCenter) / yScale;
designMatrix = polynomialDesignMatrix(u, v, fitDegree);
fitCoefficients = designMatrix \ normalizedBrightness;
fittedAtSamples = designMatrix * fitCoefficients;
fitResidual = normalizedBrightness - fittedAtSamples;
fitR2 = 1 - sum(fitResidual.^2) / ...
    max(sum((normalizedBrightness - mean(normalizedBrightness)).^2), eps);

integerX = (ceil(min(coordinateX)):floor(max(coordinateX))).';
integerY = (ceil(min(coordinateY)):floor(max(coordinateY))).';
[integerYGrid, integerXGrid] = meshgrid(integerY, integerX);
uInteger = (integerXGrid - xCenter) / xScale;
vInteger = (integerYGrid - yCenter) / yScale;
[integerDesign, ~] = polynomialDesignMatrix( ...
    uInteger(:), vInteger(:), fitDegree);
integerPrediction = reshape(integerDesign * fitCoefficients, ...
    size(integerXGrid));
[fittedIntegerMaximum, fittedIntegerIndex] = max(integerPrediction(:));
fittedIntegerX = integerXGrid(fittedIntegerIndex);
fittedIntegerY = integerYGrid(fittedIntegerIndex);
maximumIsOnBoundary = fittedIntegerX == min(integerX) || ...
    fittedIntegerX == max(integerX) || fittedIntegerY == min(integerY) || ...
    fittedIntegerY == max(integerY);

brightnessTable = table(coordinateX, coordinateY, imageCount, ...
    meanBrightness, normalizedBrightness, maxImageBrightness, ...
    folderName, brightestImage);
writetable(brightnessTable, fullfile(outputFolder, 'brightness_by_coordinate.csv'));
save(fullfile(outputFolder, 'brightness_by_coordinate.mat'), ...
    'brightnessTable', 'globalMaximum', 'globalMaximumIndex');

maxX = coordinateX(globalMaximumIndex);
maxY = coordinateY(globalMaximumIndex);
fprintf('\nGlobal maximum coordinate brightness:\n');
fprintf('  Coordinate       = (%g, %g)\n', maxX, maxY);
fprintf('  Folder           = %s\n', folderName(globalMaximumIndex));
fprintf('  Mean brightness  = %.10g\n', globalMaximum);
fprintf('  Normalized value = %.6f\n', normalizedBrightness(globalMaximumIndex));
fprintf('  Brightest image  = %s\n\n', brightestImage(globalMaximumIndex));
fprintf('Degree-%d polynomial-fit integer-coordinate maximum:\n', fitDegree);
fprintf('  Coordinate       = (%d, %d)\n', fittedIntegerX, fittedIntegerY);
fprintf('  Fitted brightness= %.8f\n', fittedIntegerMaximum);
fprintf('  Fit R^2          = %.6f\n', fitR2);
fprintf('  On scan boundary = %d\n\n', maximumIsOnBoundary);
if maximumIsOnBoundary
    warning(['The fitted integer maximum lies on the scan boundary. ', ...
        'The true peak may be outside the measured coordinate range.']);
end

%% Build the normalized 2-D coordinate grid
uniqueX = unique(coordinateX, 'sorted');
uniqueY = unique(coordinateY, 'sorted');
brightnessGrid = nan(numel(uniqueY), numel(uniqueX));
for iPoint = 1:numel(meanBrightness)
    ix = find(uniqueX == coordinateX(iPoint), 1);
    iy = find(uniqueY == coordinateY(iPoint), 1);
    brightnessGrid(iy, ix) = normalizedBrightness(iPoint);
end

%% Shared sparse horizontal ticks for coordinate figures
% Keep the fitted map and the discrete-data plots on exactly the same ticks.
% Show every second measured coordinate so that X tick labels remain horizontal.
sharedXTicks2D = uniqueY(1:2:end);
sharedYTicks2D = uniqueX(1:2:end);

% Always keep the last measured coordinate visible.
if sharedXTicks2D(end) ~= uniqueY(end)
    sharedXTicks2D(end+1) = uniqueY(end);
end
if sharedYTicks2D(end) ~= uniqueX(end)
    sharedYTicks2D(end+1) = uniqueX(end);
end

%% Curve 1: all positions in sorted scan order
figure1 = figure('Color', 'w', 'Name', 'Normalized brightness - scan order', ...
    'Units', 'centimeters', ...
    'Position', [2 2 pub.singleWidthCm pub.height2DCm]);

plot(normalizedBrightness, '-o', ...
    'LineWidth', pub.lineWidth, ...
    'MarkerSize', pub.markerSize);
hold on;
plot(globalMaximumIndex, 1, 'rp', ...
    'MarkerSize', 9, ...
    'MarkerFaceColor', 'r');

grid on;
xlabel('Coordinate index (sorted by X, then Y)', ...
    'FontName', pub.fontName, 'FontSize', pub.labelFontSize);
ylabel('Normalized Intensity', ...
    'FontName', pub.fontName, 'FontSize', pub.labelFontSize);
title('');
ylim([0, 1.08]);

text(globalMaximumIndex, 1, sprintf('  max at (%g,%g)', maxX, maxY), ...
    'Color', 'r', ...
    'FontName', pub.fontName, ...
    'FontSize', pub.tickFontSize, ...
    'FontWeight', 'bold', ...
    'VerticalAlignment', 'bottom');

set(gca, ...
    'FontName', pub.fontName, ...
    'FontSize', pub.tickFontSize, ...
    'LineWidth', pub.axisLineWidth, ...
    'TickDir', 'out', ...
    'Box', 'on');

exportgraphics(figure1, fullfile(outputFolder, ...
    'normalized_brightness_scan_order.png'), ...
    'Resolution', pub.pngResolution);
exportgraphics(figure1, fullfile(outputFolder, ...
    'normalized_brightness_scan_order.pdf'), ...
    'ContentType', 'vector', 'BackgroundColor', 'white');
savefig(figure1, fullfile(outputFolder, 'normalized_brightness_scan_order.fig'));


%% Curve 2: X-brightness curve for each Y coordinate
figure2 = figure('Color', 'w', 'Name', 'Brightness curves by coordinate', ...
    'Units', 'centimeters', ...
    'Position', [2 2 pub.mediumWidthCm pub.height2DCm]);

hold on;
colors = lines(numel(uniqueY));
for iy = 1:numel(uniqueY)
    plot(uniqueX, brightnessGrid(iy,:), '-o', ...
        'LineWidth', pub.lineWidth, ...
        'MarkerSize', pub.markerSize, ...
        'Color', colors(iy,:), ...
        'DisplayName', sprintf('Y = %g', uniqueY(iy)));
end

plot(maxX, 1, 'rp', ...
    'MarkerSize', 9, ...
    'MarkerFaceColor', 'r', ...
    'DisplayName', sprintf('Maximum (%g,%g)', maxX, maxY));

grid on;
xlabel('Coordinate X', ...
    'FontName', pub.fontName, 'FontSize', pub.labelFontSize);
ylabel('Normalized Intensity', ...
    'FontName', pub.fontName, 'FontSize', pub.labelFontSize);
title('');
ylim([0, 1.08]);

lgd = legend('Location', 'bestoutside');
set(lgd, ...
    'FontName', pub.fontName, ...
    'FontSize', pub.tickFontSize - 0.5, ...
    'Box', 'off');

set(gca, ...
    'FontName', pub.fontName, ...
    'FontSize', pub.tickFontSize, ...
    'LineWidth', pub.axisLineWidth, ...
    'TickDir', 'out', ...
    'Box', 'on');

exportgraphics(figure2, fullfile(outputFolder, ...
    'normalized_brightness_coordinate_curves.png'), ...
    'Resolution', pub.pngResolution);
exportgraphics(figure2, fullfile(outputFolder, ...
    'normalized_brightness_coordinate_curves.pdf'), ...
    'ContentType', 'vector', 'BackgroundColor', 'white');
savefig(figure2, fullfile(outputFolder, ...
    'normalized_brightness_coordinate_curves.fig'));


%% Two-dimensional heatmap in Pixel_<first>_<second> coordinate order
% Rows follow the first coordinate; columns follow the second coordinate.
% brightnessGrid is [second, first], hence .' here.
figure3 = figure('Color', 'w', 'Name', '2-D coordinate brightness heatmap', ...
    'Units', 'centimeters', ...
    'Position', [2 2 pub.singleWidthCm pub.height2DCm]);

imagesc(uniqueY, uniqueX, brightnessGrid.');

xticks(sharedXTicks2D);
yticks(sharedYTicks2D);
ax3.XTickLabelRotation = 0;
ax3.YTickLabelRotation = 0;


ax3 = gca;
% Make row coordinates increase from top to bottom (image convention)
set(ax3, 'YDir', 'reverse');
axis(ax3, 'tight');
axis(ax3, 'image');

colormap(turbo);
colorbarHandle = colorbar;
% Remove textual labels so user can add them manually later
colorbarHandle.Label.String = '';
caxis([0 1]);
colorbarHandle.Ticks = 0:0.2:1;

hold on;
% plot(maxY, maxX, 'wp', 'MarkerSize', 8, 'LineWidth', 1.2);


% Remove axis and title text for manual addition later
xlabel('', 'FontName', pub.fontName, 'FontSize', pub.labelFontSize);
ylabel('', 'FontName', pub.fontName, 'FontSize', pub.labelFontSize);
title('');

set(ax3, ...
    'FontName', pub.fontName, ...
    'FontSize', pub.tickFontSize, ...
    'LineWidth', pub.axisLineWidth, ...
    'TickDir', 'out', ...
    'Box', 'on');

set(colorbarHandle, ...
    'FontName', pub.fontName, ...
    'FontSize', pub.tickFontSize);
colorbarHandle.Label.FontName = pub.fontName;
colorbarHandle.Label.FontSize = pub.labelFontSize;
colorbarHandle.Label.Rotation = 270;
colorbarHandle.Label.HorizontalAlignment = 'center';
colorbarHandle.Label.VerticalAlignment = 'bottom';
colorbarHandle.Label.Units = 'data';

% Make the data region larger and reduce unnecessary white margin.
ax3.Position = [0.12 0.16 0.60 0.76];
colorbarHandle.Position = [0.76 0.16 0.030 0.76];

% Force horizontal tick labels after the final layout update.
drawnow;
ax3.XTickLabelRotation = 0;
ax3.YTickLabelRotation = 0;

exportgraphics(figure3, fullfile(outputFolder, ...
    'normalized_brightness_coordinate_heatmap.png'), ...
    'Resolution', pub.pngResolution);
exportgraphics(figure3, fullfile(outputFolder, ...
    'normalized_brightness_coordinate_heatmap.pdf'), ...
    'ContentType', 'vector', 'BackgroundColor', 'white');
savefig(figure3, fullfile(outputFolder, ...
    'normalized_brightness_coordinate_heatmap.fig'));


%%% Three-dimensional measured brightness distribution
% X axis: second coordinate
% Y axis: first coordinate
% Z/color: normalized brightness

[secondCoordinateGrid, firstCoordinateGrid] = meshgrid(uniqueY, uniqueX);
brightnessSurface = brightnessGrid.';

figure4 = figure( ...
    'Color', 'w', ...
    'Name', 'Measured brightness distribution', ...
    'Units', 'centimeters', ...
    'Position', [2 2 pub.singleWidthCm pub.height3DCm]);

surfaceHandle = surf( ...
    secondCoordinateGrid, ...
    firstCoordinateGrid, ...
    brightnessSurface, ...
    brightnessSurface, ...
    'FaceColor', 'interp', ...
    'EdgeColor', [0.35 0.35 0.35], ...
    'EdgeAlpha', 0.22);

colormap(turbo);

colorbarHandle = colorbar;
colorbarHandle.Label.String = 'Normalized Intensity';
colorbarHandle.Label.Rotation = 270;

caxis([0 1]);
colorbarHandle.Ticks = 0:0.2:1;
zlim([0 1]);

commonXTicks = sharedXTicks2D;
commonYTicks = sharedYTicks2D;
xticks(commonXTicks);
yticks(commonYTicks);

xlim([min(uniqueY), max(uniqueY)]);
ylim([min(uniqueX), max(uniqueX)]);

xlabel('u position', ...
    'FontName', pub.fontName, 'FontSize', pub.labelFontSize);
ylabel('v position', ...
    'FontName', pub.fontName, 'FontSize', pub.labelFontSize);
% Z axis values are normalized intensity; the colorbar provides the label.

title('');
grid on;
box on;

view(45, 30);
set(gca, 'Projection', 'perspective');
pbaspect(gca, [1 1 0.78]);

set(gca, ...
    'FontName', pub.fontName, ...
    'FontSize', pub.tickFontSize, ...
    'LineWidth', pub.axisLineWidth, ...
    'TickDir', 'out');

set(colorbarHandle, ...
    'FontName', pub.fontName, ...
    'FontSize', pub.tickFontSize);
colorbarHandle.Label.FontName = pub.fontName;
colorbarHandle.Label.FontSize = pub.labelFontSize;
colorbarHandle.Label.Rotation = 270;
colorbarHandle.Label.HorizontalAlignment = 'center';
colorbarHandle.Label.VerticalAlignment = 'bottom';
colorbarHandle.Label.Units = 'data';

% Compact publication layout: reserve room for axis labels and colorbar.
ax4 = gca;
ax4.Position = [0.08 0.15 0.62 0.74];
colorbarHandle.Position = [0.77 0.18 0.030 0.68];

material(surfaceHandle, 'dull');

exportgraphics(figure4, ...
    fullfile(outputFolder, 'measured_brightness_distribution.png'), ...
    'Resolution', pub.pngResolution);
exportgraphics(figure4, ...
    fullfile(outputFolder, 'measured_brightness_distribution.pdf'), ...
    'ContentType', 'vector', 'BackgroundColor', 'white');

savefig(figure4, ...
    fullfile(outputFolder, 'measured_brightness_distribution.fig'));


%% Polynomial fitted brightness surface

denseX = linspace(min(coordinateX), max(coordinateX), 161);
denseY = linspace(min(coordinateY), max(coordinateY), 161);

[denseYGrid, denseXGrid] = meshgrid(denseY, denseX);

uDense = (denseXGrid - xCenter) / xScale;
vDense = (denseYGrid - yCenter) / yScale;

[denseDesign, ~] = polynomialDesignMatrix( ...
    uDense(:), vDense(:), fitDegree);

densePrediction = reshape( ...
    denseDesign * fitCoefficients, ...
    size(denseXGrid));

figure5 = figure( ...
    'Color', 'w', ...
    'Name', sprintf('Degree-%d polynomial brightness fit', fitDegree), ...
    'Units', 'centimeters', ...
    'Position', [2 2 pub.singleWidthCm pub.height3DCm]);

%% Fitted surface
surfaceFit = surf( ...
    denseYGrid, ...
    denseXGrid, ...
    densePrediction, ...
    densePrediction, ...
    'EdgeColor', 'none', ...
    'FaceColor', 'interp');

hold on;

%% Original measured points
scatter3( ...
    coordinateY, ...
    coordinateX, ...
    normalizedBrightness, ...
    20, ...
    normalizedBrightness, ...
    'filled', ...
    'MarkerEdgeColor', [0.15 0.15 0.15], ...
    'LineWidth', 0.5);

colormap(turbo);

colorbarHandle = colorbar;
colorbarHandle.Label.String = 'Normalized Intensity';
colorbarHandle.Label.Rotation = 270;

caxis([0 1]);
colorbarHandle.Ticks = 0:0.2:1;
zlim([0 1]);

xticks(commonXTicks);
yticks(commonYTicks);

xlim([min(uniqueY), max(uniqueY)]);
ylim([min(uniqueX), max(uniqueX)]);

xlabel('u position', ...
    'FontName', pub.fontName, 'FontSize', pub.labelFontSize);
ylabel('v position', ...
    'FontName', pub.fontName, 'FontSize', pub.labelFontSize);
% Z axis values are normalized intensity; the colorbar provides the label.

title('');
grid on;
box on;

view(45, 32);
set(gca, 'Projection', 'perspective');

set(gca, ...
    'FontName', pub.fontName, ...
    'FontSize', pub.tickFontSize, ...
    'LineWidth', pub.axisLineWidth, ...
    'TickDir', 'out');

set(colorbarHandle, ...
    'FontName', pub.fontName, ...
    'FontSize', pub.tickFontSize);
colorbarHandle.Label.FontName = pub.fontName;
colorbarHandle.Label.FontSize = pub.labelFontSize;
colorbarHandle.Label.Rotation = 270;
colorbarHandle.Label.HorizontalAlignment = 'center';
colorbarHandle.Label.VerticalAlignment = 'bottom';
colorbarHandle.Label.Units = 'data';

ax5 = gca;
ax5.Position = [0.08 0.15 0.62 0.74];
colorbarHandle.Position = [0.77 0.18 0.030 0.68];

material(surfaceFit, 'dull');

exportgraphics(figure5, ...
    fullfile(outputFolder, 'polynomial_brightness_fit.png'), ...
    'Resolution', pub.pngResolution);
exportgraphics(figure5, ...
    fullfile(outputFolder, 'polynomial_brightness_fit.pdf'), ...
    'ContentType', 'vector', 'BackgroundColor', 'white');

savefig(figure5, ...
    fullfile(outputFolder, 'polynomial_brightness_fit.fig'));


%% Figure 6: Two-dimensional polynomial fitted brightness map
figure6 = figure( ...
    'Color', 'w', ...
    'Name', sprintf('2-D degree-%d polynomial brightness fit', fitDegree), ...
    'Units', 'centimeters', ...
    'Position', [2 2 pub.singleWidthCm pub.height2DCm]);

imagesc(denseY, denseX, densePrediction);

ax6 = gca;
% Make row coordinates increase from top to bottom (image convention)
set(ax6, 'YDir', 'reverse');
axis(ax6, 'tight');
axis(ax6, 'image');
hold on;

% Mark the fitted integer-coordinate maximum without a text annotation.
plot(fittedIntegerY, fittedIntegerX, 'wo', ...
    'MarkerSize', 5.5, ...
    'LineWidth', 1.2);


colormap(turbo);
colorbarHandle = colorbar;
% Remove textual labels so user can add them manually later
colorbarHandle.Label.String = '';
caxis([0 1]);
colorbarHandle.Ticks = 0:0.2:1;

% Use the same sparse ticks as the discrete coordinate plots.
xticks(sharedXTicks2D);
yticks(sharedYTicks2D);
ax6.XTickLabelRotation = 0;
ax6.YTickLabelRotation = 0;
xlim([min(uniqueY), max(uniqueY)]);
ylim([min(uniqueX), max(uniqueX)]);


% Remove axis text; user will add labels manually
xlabel('', 'FontName', pub.fontName, 'FontSize', pub.labelFontSize);
ylabel('', 'FontName', pub.fontName, 'FontSize', pub.labelFontSize);
title('');
box on;

set(ax6, ...
    'FontName', pub.fontName, ...
    'FontSize', pub.tickFontSize, ...
    'LineWidth', pub.axisLineWidth, ...
    'TickDir', 'out');

set(colorbarHandle, ...
    'FontName', pub.fontName, ...
    'FontSize', pub.tickFontSize);
colorbarHandle.Label.FontName = pub.fontName;
colorbarHandle.Label.FontSize = pub.labelFontSize;
colorbarHandle.Label.Rotation = 270;
colorbarHandle.Label.HorizontalAlignment = 'center';
colorbarHandle.Label.VerticalAlignment = 'bottom';
colorbarHandle.Label.Units = 'data';

% Increase the fraction of the canvas occupied by the actual data.
% These normalized positions are tuned for an 8.6 cm wide manuscript figure.
ax6.Position = [0.12 0.16 0.60 0.76];
colorbarHandle.Position = [0.76 0.16 0.030 0.76];

% Force horizontal tick labels after the final layout update.
drawnow;
ax6.XTickLabelRotation = 0;
ax6.YTickLabelRotation = 0;

exportgraphics(figure6, ...
    fullfile(outputFolder, 'polynomial_brightness_fit_2d.png'), ...
    'Resolution', pub.pngResolution);
exportgraphics(figure6, ...
    fullfile(outputFolder, 'polynomial_brightness_fit_2d.pdf'), ...
    'ContentType', 'vector', 'BackgroundColor', 'white');

savefig(figure6, ...
    fullfile(outputFolder, 'polynomial_brightness_fit_2d.fig'));


%% Figure 7: Actual scanning positions and brightest measured point

figure7 = figure( ...
    'Color', 'w', ...
    'Name', 'Actual scanning positions', ...
    'Units', 'centimeters', ...
    'Position', [2 2 pub.singleWidthCm pub.height2DCm]);

ax7 = axes(figure7);
hold(ax7, 'on');

%% 1. Plot all actual scanning positions
% coordinateX = first coordinate
% coordinateY = second coordinate
% horizontal axis -> second coordinate
% vertical axis   -> first coordinate

scatter( ...
    coordinateY, ...
    coordinateX, ...
    22, ...
    [0.45 0.12 0.12], ...
    'd', ...
    'filled', ...
    'MarkerEdgeColor', 'none');

%% 2. Brightest measured point
% If you want to display it, uncomment the block below.
%
% brightestX = coordinateX(globalMaximumIndex);
% brightestY = coordinateY(globalMaximumIndex);
% scatter( ...
%     brightestY, ...
%     brightestX, ...
%     28, ...
%     'r', ...
%     'd', ...
%     'LineWidth', 1.2, ...
%     'MarkerFaceColor', 'r', ...
%     'MarkerEdgeColor', 'none');

%% 3. Axis
xlabel('u position', ...
    'FontName', pub.fontName, 'FontSize', pub.labelFontSize);
ylabel('v position', ...
    'FontName', pub.fontName, 'FontSize', pub.labelFontSize);

xRange = max(coordinateY) - min(coordinateY);
yRange = max(coordinateX) - min(coordinateX);

xMargin = max(0.05 * xRange, 1);
yMargin = max(0.05 * yRange, 1);

xlim([ ...
    min(coordinateY) - xMargin, ...
    max(coordinateY) + xMargin]);

ylim([ ...
    min(coordinateX) - yMargin, ...
    max(coordinateX) + yMargin]);

%% 4. Use the same sparse ticks as the fitted/discrete heatmaps
ax7 = gca;
xticks(sharedXTicks2D);
yticks(sharedYTicks2D);
ax7.XTickLabelRotation = 0;
ax7.YTickLabelRotation = 0;

%% 5. Figure style
axis equal;
box on;

set(ax7, ...
    'FontName', pub.fontName, ...
    'FontSize', pub.tickFontSize, ...
    'LineWidth', pub.axisLineWidth, ...
    'TickDir', 'out');

drawnow;
ax7.XTickLabelRotation = 0;
ax7.YTickLabelRotation = 0;

title('');

%% 6. Save
exportgraphics( ...
    figure7, ...
    fullfile(outputFolder, 'actual_scan_positions.png'), ...
    'Resolution', pub.pngResolution);

exportgraphics( ...
    figure7, ...
    fullfile(outputFolder, 'actual_scan_positions.pdf'), ...
    'ContentType', 'vector', ...
    'BackgroundColor', 'white');

savefig( ...
    figure7, ...
    fullfile(outputFolder, 'actual_scan_positions.fig'));


function [design, powers] = polynomialDesignMatrix(u, v, degree)
%POLYNOMIALDESIGNMATRIX Complete 2-D polynomial basis through DEGREE.
u = u(:);
v = v(:);
termCount = (degree + 1) * (degree + 2) / 2;
design = zeros(numel(u), termCount);
powers = zeros(termCount, 2);
term = 0;
for totalDegree = 0:degree
    for powerU = totalDegree:-1:0
        term = term + 1;
        powerV = totalDegree - powerU;
        design(:, term) = u.^powerU .* v.^powerV;
        powers(term, :) = [powerU, powerV];
    end
end
end