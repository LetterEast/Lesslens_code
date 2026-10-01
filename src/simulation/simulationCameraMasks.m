function masks=simulationCameraMasks(sampleSupport,planeEdges,n)
%SIMULATIONCAMERAMASKS Mark detector pixels whose rays fall inside the image.
% These are acquisition masks, not a claim that diffraction vanishes outside
% the geometric image. Masked detector pixels are zero-filled in saved data
% and omitted from amplitude projection (not treated as measured darkness).
[r,c]=find(sampleSupport);assert(~isempty(r),'Empty sample support.');
bounds=[min(c)-.5,min(r)-.5,max(c)+.5,max(r)+.5];
masks=cell(size(planeEdges,1),1);
for k=1:numel(masks)
    b=planeEdges(k,:);
    xx=b(1)+(.5:n-.5)*(b(3)-b(1))/n;
    yy=b(2)+(.5:n-.5)*(b(4)-b(2))/n;
    [xx,yy]=meshgrid(xx,yy);
    masks{k}=xx>=bounds(1) & xx<bounds(3) & yy>=bounds(2) & yy<bounds(4);
end
end
