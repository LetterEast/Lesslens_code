function [geometry,diagnostics]=calibratePointSource(images,positions,options)
%CALIBRATEPOINTSOURCE Estimate source position from a translated detector stack.
% Refactors mainV5: projection peaks -> local centroids -> grid pitch/centre
% -> axial regression. positions are absolute offsets from FIRST detector [m].
% M,N are source offsets [pixels] from (floor(width/2)+1,floor(height/2)+1).
% Z is source-to-first-detector distance [m]. No known source is used in fitting.
if nargin<3,options=defaultCalibrationOptions();end
positions=positions(:);count=numel(images);
validateattributes(positions,{'numeric'},{'finite','real','numel',count});
assert(count>=3 && positions(1)==0 && all(diff(positions)>0), ...
    'Provide at least three increasing positions starting at zero.');
pitch=nan(count,1);centres=nan(count,2);points=cell(count,1);
for k=1:count
    image=double(images{k});if ndims(image)==3,image=image(:,:,2);end
    if k==1,imageSize=size(image);end
    assert(isequal(size(image),imageSize),'Calibration frames must have identical sizes.');
    background=imgaussfilt(image,options.backgroundSigma);
    ratio=image./max(background,eps)*mean(image(:));
    contrast=1-mat2gray(ratio);
    smooth=imgaussfilt(contrast,2);
    px=sum(smooth,1);py=sum(smooth,2);
    px=movmean(px-min(px),options.projectionSmooth);
    py=movmean(py-min(py),options.projectionSmooth);
    [vx,lx]=findProjectionPeaks(px,options.pitchPixels);
    [vy,ly]=findProjectionPeaks(py,options.pitchPixels);
    if numel(lx)<3 || numel(ly)<3,continue;end
    [~,cx]=max(vx);[~,cy]=max(vy);
    labelsX=(1:numel(lx))-cx;labelsY=(1:numel(ly))-cy;
    radius=round(options.pitchPixels/4); samples=[];
    for ix=1:numel(lx)
        for iy=1:numel(ly)
            x=round(lx(ix));y=round(ly(iy));
            if x-radius<1 || y-radius<1 || x+radius>imageSize(2) || y+radius>imageSize(1),continue;end
            xx=x-radius:x+radius;yy=y-radius:y+radius;
            patch=contrast(yy,xx);patch=patch-min(patch(:));weight=sum(patch(:));
            if weight<=eps,continue;end
            [gx,gy]=meshgrid(xx,yy);
            samples(end+1,:)=[labelsX(ix),labelsY(iy), ...
                sum(gx(:).*patch(:))/weight,sum(gy(:).*patch(:))/weight]; %#ok<AGROW>
        end
    end
    if size(samples,1)<6,continue;end
    fitX=polyfit(samples(:,1),samples(:,3),1);fitY=polyfit(samples(:,2),samples(:,4),1);
    pitch(k)=mean([fitX(1),fitY(1)]);centres(k,:)=[fitX(2),fitY(2)];points{k}=samples;
    if options.showFigures
        figure(731);imshow(image,[]);hold on;plot(samples(:,3),samples(:,4),'r.');
        plot(centres(k,1),centres(k,2),'g+');hold off;title(sprintf('Calibration frame %d/%d',k,count));drawnow;
    end
end
valid=isfinite(pitch)&pitch>0;
assert(nnz(valid)>=3,'Too few valid calibration frames; check grid pitch and image contrast.');
% Keep ORIGINAL positions when frames fail; never renumber accepted frames.
p=polyfit(positions(valid),pitch(valid),1);
assert(p(1)>0 && p(2)>0,'Grid pitch must increase with positive detector displacement.');
magnification=pitch(valid)/p(2);
fitX=polyfit(magnification,centres(valid,1),1);
fitY=polyfit(magnification,centres(valid,2),1);
geometry=struct('M',fitX(2)-(floor(imageSize(2)/2)+1), ...
    'N',fitY(2)-(floor(imageSize(1)/2)+1),'Z',p(2)/p(1));
diagnostics=struct('positions',positions,'validFrames',valid,'pitchPixels',pitch, ...
    'centresPixels',centres,'points',{points},'pitchFit',p, ...
    'centreFitX',fitX,'centreFitY',fitY,'options',options, ...
    'pitchResidualPixels',pitch(valid)-polyval(p,positions(valid)));
end
