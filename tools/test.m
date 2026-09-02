clc;
close all;
clear;

%% Configuration
dataRoot = "\\192.168.2.166\d\lesslens\2026.9.1\NONE_B";
imageSubfolder = fullfile("pixel_on", "exposure_-1");
outputFolder = fullfile(fileparts(mfilename('fullpath')), ...
    'ResultFolder', 'BrightnessAnalysis_20260901');
validExtensions = [".png", ".jpg", ".jpeg", ".bmp", ".tif", ".tiff"];

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

%% Fit a 2-D quadratic surface and search only integer coordinates
% Center/scale coordinates to improve numerical conditioning. The fitted
% model is z = b0+b1*u+b2*v+b3*u^2+b4*u*v+b5*v^2.
xCenter = mean(coordinateX);
yCenter = mean(coordinateY);
xScale = max(std(coordinateX), 1);
yScale = max(std(coordinateY), 1);
u = (coordinateX - xCenter) / xScale;
v = (coordinateY - yCenter) / yScale;
designMatrix = [ones(size(u)), u, v, u.^2, u.*v, v.^2];
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
integerDesign = [ones(numel(uInteger),1), uInteger(:), vInteger(:), ...
    uInteger(:).^2, uInteger(:).*vInteger(:), vInteger(:).^2];
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
fprintf('Quadratic-fit integer-coordinate maximum:\n');
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

%% Curve 1: all positions in sorted scan order
figure1 = figure('Color', 'w', 'Name', 'Normalized brightness - scan order');
plot(normalizedBrightness, '-o', 'LineWidth', 1.4, 'MarkerSize', 4);
hold on;
plot(globalMaximumIndex, 1, 'rp', 'MarkerSize', 13, 'MarkerFaceColor', 'r');
grid on;
xlabel('Coordinate index (sorted by X, then Y)');
ylabel('Normalized mean brightness');
title('Normalized brightness of all coordinate positions');
ylim([0, 1.08]);
text(globalMaximumIndex, 1, sprintf('  max at (%g,%g)', maxX, maxY), ...
    'Color', 'r', 'FontWeight', 'bold', 'VerticalAlignment', 'bottom');
exportgraphics(figure1, fullfile(outputFolder, ...
    'normalized_brightness_scan_order.png'), 'Resolution', 200);
savefig(figure1, fullfile(outputFolder, 'normalized_brightness_scan_order.fig'));

%% Curve 2: X-brightness curve for each Y coordinate
figure2 = figure('Color', 'w', 'Name', 'Brightness curves by coordinate');
hold on;
colors = lines(numel(uniqueY));
for iy = 1:numel(uniqueY)
    plot(uniqueX, brightnessGrid(iy,:), '-o', 'LineWidth', 1.5, ...
        'MarkerSize', 5, 'Color', colors(iy,:), ...
        'DisplayName', sprintf('Y = %g', uniqueY(iy)));
end
plot(maxX, 1, 'rp', 'MarkerSize', 14, 'MarkerFaceColor', 'r', ...
    'DisplayName', sprintf('Maximum (%g,%g)', maxX, maxY));
grid on;
xlabel('Coordinate X');
ylabel('Normalized mean brightness');
title('Normalized brightness versus coordinate X');
ylim([0, 1.08]);
legend('Location', 'bestoutside');
exportgraphics(figure2, fullfile(outputFolder, ...
    'normalized_brightness_coordinate_curves.png'), 'Resolution', 200);
savefig(figure2, fullfile(outputFolder, ...
    'normalized_brightness_coordinate_curves.fig'));

%% Two-dimensional heatmap in Pixel_<first>_<second> coordinate order
% Rows follow the first coordinate (194 -> 204); columns follow the second
% coordinate (305 -> 310). brightnessGrid is [second, first], hence .' here.
figure3 = figure('Color', 'w', 'Name', '2-D coordinate brightness heatmap');
imagesc(uniqueY, uniqueX, brightnessGrid.');
set(gca, 'YDir', 'normal');
axis tight;
xticks(uniqueY);
yticks(uniqueX);
colormap(turbo);
colorbarHandle = colorbar;
colorbarHandle.Label.String = 'Normalized mean brightness';
caxis([0 1]);
hold on;
plot(maxY, maxX, 'wp', 'MarkerSize', 15, 'LineWidth', 1.5, ...
    'MarkerFaceColor', 'r');
text(maxY, maxX, sprintf('  max %.4g', globalMaximum), ...
    'Color', 'w', 'FontWeight', 'bold', 'VerticalAlignment', 'bottom');
xlabel('Second coordinate');
ylabel('First coordinate');
title('Normalized brightness heatmap in Pixel_{first,second} order');
exportgraphics(figure3, fullfile(outputFolder, ...
    'normalized_brightness_coordinate_heatmap.png'), 'Resolution', 200);
savefig(figure3, fullfile(outputFolder, ...
    'normalized_brightness_coordinate_heatmap.fig'));

%% Three-dimensional coordinate-brightness heatmap
% X axis: second coordinate; Y axis: first coordinate; Z/color: brightness.
[secondCoordinateGrid, firstCoordinateGrid] = meshgrid(uniqueY, uniqueX);
brightnessSurface = brightnessGrid.';
figure4 = figure('Color', 'w', 'Name', '3-D coordinate brightness heatmap');
surfaceHandle = surf(secondCoordinateGrid, firstCoordinateGrid, ...
    brightnessSurface, brightnessSurface, ...
    'FaceColor', 'interp', 'EdgeColor', [0.25 0.25 0.25], ...
    'EdgeAlpha', 0.35);
hold on;
plot3(maxY, maxX, 1, 'rp', 'MarkerSize', 15, ...
    'MarkerFaceColor', 'r', 'LineWidth', 1.5);
text(maxY, maxX, 1, sprintf('  max (%g,%g) = 1', maxX, maxY), ...
    'Color', 'r', 'FontWeight', 'bold', 'VerticalAlignment', 'bottom');
colormap(turbo);
colorbarHandle = colorbar;
colorbarHandle.Label.String = 'Normalized mean brightness';
caxis([0 1]);
zlim([0 1.08]);
xticks(uniqueY);
yticks(uniqueX);
xlabel('Second coordinate');
ylabel('First coordinate');
zlabel('Normalized mean brightness');
title('3-D normalized brightness heatmap');
grid on;
box on;
view(45, 32);
set(gca, 'Projection', 'perspective');
material(surfaceHandle, 'dull');
exportgraphics(figure4, fullfile(outputFolder, ...
    'normalized_brightness_3d_heatmap.png'), 'Resolution', 240);
savefig(figure4, fullfile(outputFolder, ...
    'normalized_brightness_3d_heatmap.fig'));

%% Fitted surface and integer-constrained maximum
denseX = linspace(min(coordinateX), max(coordinateX), 161);
denseY = linspace(min(coordinateY), max(coordinateY), 161);
[denseYGrid, denseXGrid] = meshgrid(denseY, denseX);
uDense = (denseXGrid - xCenter) / xScale;
vDense = (denseYGrid - yCenter) / yScale;
denseDesign = [ones(numel(uDense),1), uDense(:), vDense(:), ...
    uDense(:).^2, uDense(:).*vDense(:), vDense(:).^2];
densePrediction = reshape(denseDesign * fitCoefficients, size(denseXGrid));

figure5 = figure('Color', 'w', 'Name', 'Quadratic brightness fit');
surf(denseYGrid, denseXGrid, densePrediction, densePrediction, ...
    'EdgeColor', 'none', 'FaceColor', 'interp');
hold on;
scatter3(coordinateY, coordinateX, normalizedBrightness, 28, ...
    normalizedBrightness, 'filled', 'MarkerEdgeColor', 'k');
plot3(fittedIntegerY, fittedIntegerX, fittedIntegerMaximum, 'rp', ...
    'MarkerSize', 16, 'MarkerFaceColor', 'r', 'LineWidth', 1.5);
text(fittedIntegerY, fittedIntegerX, fittedIntegerMaximum, ...
    sprintf('  integer max (%d,%d)', fittedIntegerX, fittedIntegerY), ...
    'Color', 'r', 'FontWeight', 'bold');
colormap(turbo);
colorbar;
xlabel('Second coordinate');
ylabel('First coordinate');
zlabel('Fitted normalized brightness');
title(sprintf('Quadratic surface fit, R^2 = %.4f', fitR2));
grid on;
view(45, 32);
exportgraphics(figure5, fullfile(outputFolder, ...
    'quadratic_fit_integer_maximum.png'), 'Resolution', 240);
savefig(figure5, fullfile(outputFolder, ...
    'quadratic_fit_integer_maximum.fig'));

fitResult = struct('coefficients', fitCoefficients, 'R2', fitR2, ...
    'integerMaximumX', fittedIntegerX, ...
    'integerMaximumY', fittedIntegerY, ...
    'integerMaximumBrightness', fittedIntegerMaximum, ...
    'maximumIsOnBoundary', maximumIsOnBoundary);
save(fullfile(outputFolder, 'quadratic_fit_result.mat'), 'fitResult');

fprintf('Processed %d coordinate folders. Results saved to:\n  %s\n', ...
    height(brightnessTable), outputFolder);
