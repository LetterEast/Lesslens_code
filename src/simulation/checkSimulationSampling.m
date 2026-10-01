function diagnostics=checkSimulationSampling(o)
%CHECKSIMULATIONSAMPLING Guard the DIRECT complex-field ASM implementation.
% This is a numerical carrier-sampling limit, not a universal camera-angle
% limit: an envelope/shifted-spectrum solver can represent larger tilts.
n=o.imageSize;p=o.padding;side=n+2*p;dx=o.pixelSize;L=o.sourceToSample;
validateattributes(L,{'numeric'},{'scalar','real','finite','positive'});
validateattributes(dx,{'numeric'},{'scalar','real','finite','positive'});
validateattributes(o.wavelength,{'numeric'},{'scalar','real','finite','positive'});
validateattributes(o.sourceOffsetPixels,{'numeric'},{'vector','numel',2,'finite','real'});
[x,y]=meshgrid(((1:side)-floor(side/2)-1)*dx);
x=x-o.sourceOffsetPixels(1)*dx;y=y-o.sourceOffsetPixels(2)*dx;
r=sqrt(x.^2+y.^2+L^2);
cyclesX=abs(x./r)*dx/o.wavelength;cyclesY=abs(y./r)*dx/o.wavelength;
active=true(side);
if isfield(o,'sampleImageSize') && strcmp(o.sample.type,'image')
    % Finite images are identically zero in the padding. Only their support
    % and the independently generated calibration target carry a wavefield.
    [gx,gy]=meshgrid((1:side)-floor(side/2)-1);
    sz=o.sampleImageSize;shift=o.sample.offsetPixels;
    active=abs(gx-shift(1))<=sz(2)/2+1 & abs(gy-shift(2))<=sz(1)/2+1;
    calHalf=o.calibration.imageSize/2;
    active=active | (abs(gx)<=calHalf & abs(gy)<=calHalf);
end
diagnostics.maxCarrierCyclesPerPixel=[max(cyclesX(active)),max(cyclesY(active))];
diagnostics.centreTiltDegrees=atan(norm(o.sourceOffsetPixels)*dx/L)*180/pi;
diagnostics.nyquistCyclesPerPixel=.5;
if max(diagnostics.maxCarrierCyclesPerPixel)>=.5
    error('Lensless:SimulationAliasing', ...
        ['Direct complex-field sampling fails: max carrier [%.3f, %.3f] cycles/pixel (limit 0.5), centre tilt %.2f deg. ' ...
        'Reduce sourceOffsetPixels, increase sourceToSample, or use a finer computation grid with a proper camera model. ' ...
        'Padding alone does not fix carrier aliasing. Restore a consistent set of parameters in main/demo_sim.m.'], ...
        diagnostics.maxCarrierCyclesPerPixel(1),diagnostics.maxCarrierCyclesPerPixel(2),diagnostics.centreTiltDegrees);
end
% A necessary (not sufficient) check: tracked calibration centre must remain
% on the finite detector; otherwise another grid dot can be mistaken for it.
calCentres=-(o.calibrationDistances(:)/L)*o.sourceOffsetPixels(:).';
diagnostics.calibrationCentreOffsets=calCentres;
calSize=n;
if isfield(o,'calibration') && isfield(o.calibration,'imageSize'),calSize=o.calibration.imageSize;end
if any(abs(calCentres(:))>calSize/2-8)
    error('Lensless:CalibrationOutsideSensor', ...
        'Calibration centre leaves the camera. Reduce source offset/calibration distances or enlarge the sensor.');
end
% Check footprint bounds before sampleGeometrySupport clips them to the grid.
source=[floor(side/2)+1,floor(side/2)+1]+o.sourceOffsetPixels;
sensor=[p+.5,p+.5,p+n+.5,p+n+.5];
edges=zeros(numel(o.sampleDistances),4);
for k=1:numel(o.sampleDistances)
    scale=L/(L+o.sampleDistances(k));
    edges(k,:)=[source source]+scale*(sensor-[source source]);
end
diagnostics.objectPlaneEdges=edges;
if any(edges(:,1:2)<.5,'all') || any(edges(:,3:4)>side+.5,'all')
    error('Lensless:SimulationFOVClipped', ...
        'A sample-plane camera footprint leaves the grid. Increase padding and recheck carrier sampling.');
end
end
