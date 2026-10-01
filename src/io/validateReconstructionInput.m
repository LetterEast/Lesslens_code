function validateReconstructionInput(inputData)
%VALIDATERECONSTRUCTIONINPUT Check the public prepared-measurement contract.
% Acquisition, file naming and raw-image registration belong to caller adapters.
required = {'images','geometry','distanceSteps','wavelength','pixelSize'};
if ~isstruct(inputData) || ~isscalar(inputData) || ~all(isfield(inputData, required))
    error('Lensless:InvalidInput', 'inputData must contain: %s.', strjoin(required, ', '));
end
images = inputData.images;
if ~iscell(images) || numel(images) < 2
    error('Lensless:InvalidInput', 'At least two prepared intensity images are required.');
end
imageSize = size(images{1});
for index = 1:numel(images)
    validateattributes(images{index}, {'single','double'}, ...
        {'real','finite','nonnegative','nonempty','2d'}, mfilename, 'images');
    if ~isequal(size(images{index}), imageSize)
        error('Lensless:InvalidInput', 'All prepared images must have the same size.');
    end
end
validateattributes(inputData.distanceSteps, {'numeric'}, ...
    {'real','finite','vector','numel',numel(images)});
validateattributes(inputData.wavelength, {'numeric'}, {'real','finite','scalar','positive'});
validateattributes(inputData.pixelSize, {'numeric'}, {'real','finite','scalar','positive'});
g = inputData.geometry;
fields = {'M','N','Z','orig_M','orig_N','padSize','orig_size','ValidMaskHard'};
if ~isstruct(g) || ~isscalar(g) || ~all(isfield(g, fields))
    error('Lensless:InvalidGeometry', 'geometry must contain: %s.', strjoin(fields, ', '));
end
for name = {'M','N','Z','orig_M','orig_N'}
    validateattributes(g.(name{1}), {'numeric'}, {'real','finite','scalar'});
end
validateattributes(g.Z, {'numeric'}, {'positive'});
validateattributes(g.padSize, {'numeric'}, {'real','finite','integer','nonnegative','numel',2});
validateattributes(g.orig_size, {'numeric'}, {'real','finite','integer','positive','numel',2});
if ~isequal(double(g.orig_size(:).'+2*g.padSize(:).'), double(imageSize))
    error('Lensless:InvalidGeometry', 'orig_size + 2*padSize must match the prepared image size.');
end
if ~iscell(g.ValidMaskHard) || numel(g.ValidMaskHard) ~= numel(images)
    error('Lensless:InvalidGeometry', 'One hard support mask per image is required.');
end
for index = 1:numel(images)
    mask = g.ValidMaskHard{index};
    if ~isequal(size(mask), imageSize) || ~any(mask(:)) || ...
            ~all(mask(:)==0 | mask(:)==1)
        error('Lensless:InvalidGeometry', 'Each mask must be binary, nonempty and match its image.');
    end
end
end
