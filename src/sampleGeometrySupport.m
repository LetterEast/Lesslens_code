function support = sampleGeometrySupport(imageSize, geometry, zPositions, focusDistance)
%SAMPLEGEOMETRYSUPPORT Geometric sensor footprints in the sample-plane grid.
% Same point-source, parallel-plane mapping as test_geometry_crop. Coordinates
% include registration shear; do not add the original source offset twice.
required = {'Z','M','N','orig_M','orig_N','padSize','orig_size'};
assert(all(isfield(geometry, required)), 'Missing geometry for sample-plane cropping.');
g = geometry;
zPositions = double(zPositions(:));
sourceToObject = double(g.Z)-double(focusDistance);
assert(double(g.Z)>0 && sourceToObject>0 && all(double(g.Z)+zPositions>0), ...
    'Source-to-object and source-to-detector distances must be positive.');
height = imageSize(1); width = imageSize(2);
source = [floor(width/2)+1, floor(height/2)+1]+double([g.M,g.N]);
offsets = double([g.orig_M,g.orig_N]);
pad = double(g.padSize); originalSize = double(g.orig_size);
edges = [pad(2)+0.5,pad(1)+0.5, ...
    pad(2)+originalSize(2)+0.5,pad(1)+originalSize(1)+0.5];
support.version = 1;
support.sourceToObject = sourceToObject;
support.magnifications = (double(g.Z)+zPositions)/sourceToObject;
support.scales = 1./support.magnifications;
support.planeEdges = zeros(numel(zPositions),4);
support.coverage = zeros(height,width,'uint16');
support.firstMask = false(height,width);
for index = 1:numel(zPositions)
    shift = zPositions(index)*offsets/double(g.Z);
    b = [source,source]+support.scales(index)*(edges+[shift,shift]-[source,source]);
    support.planeEdges(index,:) = b;
    columns = max(1,ceil(b(1))):min(width,floor(b(3)));
    rows = max(1,ceil(b(2))):min(height,floor(b(4)));
    support.coverage(rows,columns) = support.coverage(rows,columns)+1;
    if index==1, support.firstMask(rows,columns) = true; end
end
support.mask = support.coverage>0;
[rows,columns] = find(support.mask);
assert(~isempty(rows), 'Geometric sample footprint is empty.');
assert(any(support.firstMask(:)), 'First-camera sample footprint is empty.');
support.bounds = [min(columns),min(rows),max(columns),max(rows)];
end
