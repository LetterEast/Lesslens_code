function [truth,sampleSupport]=createSimulationSample(n,p,sample)
%CREATESIMULATIONSAMPLE Built-in analytic objects or user amplitude/phase images.
% Images are converted to grayscale [0,1], aspect-preserving resized into
% sample.sizePixels along the longest side, independent of detector width.
% amplitude = 1 - absorptionScale*(1-gray), clipped at zero.
% phase = phaseScale*gray radians. No per-image contrast stretching is used.
% Built-ins are generated here and need no downloaded assets.
validateattributes(sample.phaseScale,{'numeric'},{'scalar','finite','real','nonnegative'});
validateattributes(sample.absorptionScale,{'numeric'},{'scalar','finite','real','nonnegative'});
side=n+2*p;[x,y]=meshgrid((1:side)-floor(side/2)-1);
sampleSupport=true(side);
envelope=exp(-((x/(.29*n)).^8+(y/(.29*n)).^8));
phase=zeros(side);absorption=zeros(side);
switch validatestring(sample.type,{'mixed','resolution','cells','image','widefield'})
    case 'widefield'
        % Repeated absorbing dots and phase bumps extend beyond the first FOV.
        for cx=-.55*n:.11*n:.55*n
            for cy=-.55*n:.11*n:.55*n
                r2=(x-cx).^2+(y-cy).^2;
                absorption=absorption+.55*exp(-r2/(2*(.022*n)^2));
                phase=phase+.6*exp(-r2/(2*(.030*n)^2));
            end
        end
        absorption=min(.95,sample.absorptionScale*absorption);
    case 'mixed'
        phase=(.9*exp(-((x+18).^2+(y+12).^2)/250)+ ...
            .55*exp(-((x-20).^2+(y-16).^2)/160)).*envelope;
        disk=.3*exp(-((x+20).^2+(y-18).^2)/90);
        bars=(abs(y+20)<10)&(x>0)&(x<40)&(mod(floor(x),8)<4);
        absorption=min(.65,(disk+.4*double(bars)).*envelope);
        absorption=min(.95,sample.absorptionScale*absorption);
    case 'resolution'
        % Three stripes in each orientation, with widths 1..6 pixels.
        for k=1:6
            w=k;cx=(mod(k-1,3)-1)*.22*n;cy=(floor((k-1)/3)-.5)*.3*n;
            for j=-1:1
                if mod(k,2)==1
                    dx=x-cx-2*j*w;
                    stripe=dx>=-w/2 & dx<w/2 & abs(y-cy)<.065*n;
                else
                    dy=y-cy-2*j*w;
                    stripe=dy>=-w/2 & dy<w/2 & abs(x-cx)<.065*n;
                end
                absorption(stripe)=.8;
            end
        end
        absorption=sample.absorptionScale*absorption;
    case 'cells'
        centres=[-.17 -.12 .075;.12 -.1 .09;-.03 .15 .10];
        for k=1:size(centres,1)
            r=hypot(x/n-centres(k,1),y/n-centres(k,2))/centres(k,3);
            phase=phase+.8*sqrt(max(0,1-r.^2));
        end
    case 'image'
        sampleSupport=false(side);
        assert(~isempty(sample.amplitudeFile)||~isempty(sample.phaseFile), ...
            'Image mode requires amplitudeFile or phaseFile.');
        if ~isempty(sample.amplitudeFile)
            [gray,rows,cols]=readPatch(sample.amplitudeFile,side,sample.sizePixels,sample.offsetPixels);
            absorption(rows,cols)=sample.absorptionScale*(1-gray);
            sampleSupport(rows,cols)=true;
        end
        if ~isempty(sample.phaseFile)
            [gray,rows,cols]=readPatch(sample.phaseFile,side,sample.sizePixels,sample.offsetPixels);
            phase(rows,cols)=gray;
            if isempty(sample.amplitudeFile),sampleSupport(rows,cols)=true;end
        end
end
truth=max(0,1-absorption).*exp(1i*sample.phaseScale*phase);
truth(~sampleSupport)=0; % Image boundary is finite; no white/periodic extension.
end

function [gray,rows,cols]=readPatch(file,side,sizePixels,offset)
[raw,map,alpha]=imread(file);
if ~isempty(map),gray=ind2gray(raw,map);
elseif size(raw,3)==3,gray=im2double(rgb2gray(raw));
else,gray=im2double(raw);end
assert(ismatrix(gray)&&all(isfinite(gray(:)))&&all(gray(:)>=0 & gray(:)<=1), ...
    'Sample image must have finite grayscale values in [0,1].');
assert(isempty(alpha)||all(im2double(alpha(:))==1), ...
    'Use an opaque sample image; flatten transparency before importing.');
validateattributes(sizePixels,{'numeric'},{'scalar','finite','integer','positive'});
assert(sizePixels<side-10,'Lensless:SampleOutsideGrid', ...
    'Sample longest side is %d pixels, but the grid is %d pixels. Increase padding or reduce sample.sizePixels (leave >10 pixels total margin).',sizePixels,side);
scale=sizePixels/max(size(gray));
gray=imresize(gray,max(1,round(size(gray)*scale)),'bilinear');
validateattributes(offset,{'numeric'},{'real','finite','integer','numel',2});
rows=floor((side-size(gray,1))/2)+(1:size(gray,1))+offset(2);
cols=floor((side-size(gray,2))/2)+(1:size(gray,2))+offset(1);
assert(min(rows)>5 && max(rows)<side-5 && min(cols)>5 && max(cols)<side-5, ...
    'Lensless:SampleOutsideGrid', ...
    'Shifted sample exceeds the %d-pixel grid. Increase padding or reduce sample.sizePixels / sample.offsetPixels.',side);
end
