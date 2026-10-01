function o=resolveSimulationCamera(o)
%RESOLVESIMULATIONCAMERA Keep the square camera within the resized image.
validateattributes(o.imageSize,{'numeric'},{'scalar','integer','positive','finite'});
o.requestedImageSize=o.imageSize;
if ~strcmp(o.sample.type,'image'),return;end
file=o.sample.amplitudeFile;
if isempty(file),file=o.sample.phaseFile;end
assert(~isempty(file),'Image mode requires amplitudeFile or phaseFile.');
info=imfinfo(file);shape=double([info(1).Height,info(1).Width]);
validateattributes(o.sample.sizePixels,{'numeric'},{'scalar','integer','positive','finite'});
o.sampleImageSize=max(1,round(shape*o.sample.sizePixels/max(shape)));
% The current camera is square; for rectangular images use the shorter side.
o.imageSize=min(o.imageSize,min(o.sampleImageSize));
if o.imageSize<o.requestedImageSize
    fprintf('Camera limited from %d to %d pixels; resized sample is %d x %d.\n', ...
        o.requestedImageSize,o.imageSize,o.sampleImageSize(1),o.sampleImageSize(2));
end
end
