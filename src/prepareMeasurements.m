function [images, geometry] = prepareMeasurements( ...
        imageFolder, geometry, distanceSteps)
%PREPAREMEASUREMENTS Load and register reconstruction-ready intensity data.

images = loadImages(imageFolder);
count = numel(images);
if numel(distanceSteps) ~= count
    error('Distance steps must contain one value per image.');
end
assert(all(cellfun(@(x) isequal(size(x), size(images{1})), images)), ...
    'All input images must have the same size.');

[height, width] = size(images{1});
positions = cumsum(distanceSteps(:).');
shiftX = positions * geometry.M / geometry.Z;
shiftY = positions * geometry.N / geometry.Z;
padX = ceil(max(abs(shiftX))) + 200;
padY = ceil(max(abs(shiftY))) + 200;
canvasSize = [height + 2*padY, width + 2*padX];

registered = cell(1, count);
hardMasks = cell(1, count);
softMasks = cell(1, count);
for index = 1:count
    padded = padarray(images{index}, [padY, padX], 'replicate', 'both');
    registered{index} = imtranslate(padded, [shiftX(index), shiftY(index)], ...
        'bicubic', 'OutputView', 'same', 'FillValues', 0);
    registered{index} = max(registered{index}, 0);

    firstColumn = padX + 1 + round(shiftX(index));
    firstRow = padY + 1 + round(shiftY(index));
    bounds = [firstColumn, firstRow, ...
        firstColumn + width - 1, firstRow + height - 1];
    hardMasks{index} = rectangleMask(canvasSize, bounds);
    softMasks{index} = softenMask(hardMasks{index}, 100);
end

geometry.orig_M = geometry.M;
geometry.orig_N = geometry.N;
geometry.orig_size = [height, width];
geometry.padSize = [padY, padX];
geometry.ValidMask = softMasks;
geometry.ValidMaskHard = hardMasks;
geometry.CoverageCount = sum(cat(3, hardMasks{:}), 3);
geometry.UnionMask = geometry.CoverageCount >= 1;
geometry.TrustedMask = geometry.CoverageCount >= 2;
geometry.M = 0;
geometry.N = 0;
images = registered;
end

function images = loadImages(folder)
if ~isfolder(folder), error('Image folder does not exist: %s', folder); end
entries = dir(folder);
entries = entries(~[entries.isdir]);
extensions = {'.png', '.jpg', '.jpeg', '.bmp', '.tif', '.tiff'};
keep = false(size(entries));
for index = 1:numel(entries)
    [~, ~, extension] = fileparts(entries(index).name);
    keep(index) = any(strcmpi(extension, extensions));
end
entries = entries(keep);
if isempty(entries), error('No supported images found in: %s', folder); end
[~, order] = sort(lower(string({entries.name})));
entries = entries(order);

images = cell(numel(entries), 1);
for index = 1:numel(entries)
    image = im2double(imread(fullfile(entries(index).folder, entries(index).name)));
    if size(image, 3) == 3, image = rgb2gray(image); end
    images{index} = image;
end
end

function mask = rectangleMask(canvasSize, bounds)
mask = false(canvasSize);
x1 = max(1, bounds(1));
y1 = max(1, bounds(2));
x2 = min(canvasSize(2), bounds(3));
y2 = min(canvasSize(1), bounds(4));
if x1 <= x2 && y1 <= y2, mask(y1:y2, x1:x2) = true; end
end

function mask = softenMask(hardMask, transitionWidth)
% Separable cosine edge taper inside the measured rectangle.
[rows, columns] = find(hardMask);
mask = zeros(size(hardMask));
if isempty(rows), return; end
y1 = min(rows); y2 = max(rows); x1 = min(columns); x2 = max(columns);
height = y2 - y1 + 1; width = x2 - x1 + 1;
edgeY = min([transitionWidth, floor(height/2)]);
edgeX = min([transitionWidth, floor(width/2)]);
windowY = ones(height, 1); windowX = ones(1, width);
if edgeY > 0
    ramp = 0.5 * (1 - cos(pi * (0:edgeY-1) / edgeY));
    windowY(1:edgeY) = ramp; windowY(end-edgeY+1:end) = fliplr(ramp);
end
if edgeX > 0
    ramp = 0.5 * (1 - cos(pi * (0:edgeX-1) / edgeX));
    windowX(1:edgeX) = ramp; windowX(end-edgeX+1:end) = fliplr(ramp);
end
mask(y1:y2, x1:x2) = windowY * windowX;
end
