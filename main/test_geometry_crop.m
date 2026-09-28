function outputFolder = test_geometry_crop(iterationFolder)
%TEST_GEOMETRY_CROP Project detector footprints onto the sample plane.
% Uses the geometry of the saved reconstruction, not current input settings.
% Registered coordinates are sheared by x'=x+z*orig_M/Z (and y likewise),
% so the source is at geometry.M/N (normally zero) in this coordinate frame.
% This is a ray-geometric footprint, not a diffraction reliability boundary.
projectRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(projectRoot, 'src'));
if nargin < 1 || isempty(iterationFolder)
    files = dir(fullfile(projectRoot, 'ResultFolder', 'APRW_*', ...
        'Iter_*', 'results', 'reconstruction.mat'));
    if isempty(files), error('No saved reconstruction found.'); end
    [~, index] = max([files.datenum]);
    iterationFolder = fileparts(files(index).folder);
end
r = load(fullfile(iterationFolder, 'results', 'reconstruction.mat'), ...
    'field', 'focusDistance', 'bounds');
metadata = load(fullfile(iterationFolder, 'diagnostics', 'meta.mat'), ...
    'mainPara', 'zPositions');
p = metadata.mainPara;
g = p.MNZ_result;
zPositions = double(metadata.zPositions(:));
focusDistance = double(r.focusDistance);
support = sampleGeometrySupport(size(r.field), g, zPositions, focusDistance);
sourceToObject = support.sourceToObject;
magnifications = support.magnifications;
scales = support.scales;
planeEdges = support.planeEdges;
sampleCoverage = support.coverage;
sampleMask = support.mask;
geometryBounds = support.bounds;
oldBounds = double(r.bounds);

% Keep spherical illumination during propagation. Its subsequent phase-only
% removal does not change amplitude, so no sample-phase operation is needed.
objectAmplitude = abs(propagate(r.field, p.PixelSize, p.WaveLength, -focusDistance));
oldRows = oldBounds(2):oldBounds(4);
oldColumns = oldBounds(1):oldBounds(3);
baseline = objectAmplitude(oldRows,oldColumns);
displayLimits = prctile(double(baseline(:)), [1,99]);
displayImage = min(max((double(objectAmplitude)-displayLimits(1)) / ...
    max(diff(displayLimits),eps),0),1);
newRows = geometryBounds(2):geometryBounds(4);
newColumns = geometryBounds(1):geometryBounds(3);
outputFolder = fullfile(iterationFolder, 'geometry_crop_test', ...
    char(datetime('now','Format','yyyyMMdd_HHmmss_SSS')));
mkdir(outputFolder);
imwrite(displayImage(oldRows,oldColumns), ...
    fullfile(outputFolder,'baseline_current_crop_noTV.png'));
writeMaskedPNG(displayImage(newRows,newColumns), ...
    fullfile(outputFolder,'geometry_union_crop_noTV.png'), ...
    sampleMask(newRows,newColumns));
imwrite(uint8(sampleMask(newRows,newColumns))*255, ...
    fullfile(outputFolder,'geometry_union_mask.png'));
save(fullfile(outputFolder,'geometry.mat'), 'magnifications', 'scales', ...
    'planeEdges', 'sampleCoverage', 'geometryBounds', 'oldBounds', ...
    'focusDistance', 'sourceToObject', 'zPositions', 'displayLimits', ...
    'iterationFolder');
assert(all(sampleCoverage(:)<=numel(zPositions)));
assert(all(isfinite(planeEdges(:))));
fprintf('Magnifications: '); fprintf('%.6f ',magnifications); fprintf('\n');
fprintf('Old bounds [left top right bottom]: %g %g %g %g\n',oldBounds);
fprintf('Geometric bounds: %g %g %g %g\n',geometryBounds);
fprintf('Separate images saved to: %s\n',outputFolder);
end
