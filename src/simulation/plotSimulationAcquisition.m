function plotSimulationAcquisition(result,folder)
%PLOTSIMULATIONACQUISITION Large object -> finite camera aperture at each z.
% Top panels show geometric footprints, not intensity patches cropped from
% the object: detector images below are obtained AFTER wave propagation.
o=result.options;side=size(result.truth,1);count=numel(result.images);
coord=((1:side)-floor(side/2)-1)*o.pixelSize*1e6;
edges=(result.fov.truth.planeEdges-floor(side/2)-1)*o.pixelSize*1e6;
visible='off';if o.showFigures,visible='on';end
selected=unique(round(linspace(1,count,min(6,count))));columns=numel(selected);
sampleMask=result.sampleSupport;[sr,sc]=find(sampleMask|result.fov.truth.mask);
viewX=coord([max(1,min(sc)-10),min(side,max(sc)+10)]);
viewY=coord([max(1,min(sr)-10),min(side,max(sr)+10)]);
f=figure('Visible',visible,'Name','仿真采集：大样品与有限相机','Position',[70 70 1500 950]);
for j=1:columns
    k=selected(j);
    subplot(3,columns,j);objectPanel(k);title(sprintf('Frame %d | z=%.2f mm',k,o.sampleDistances(k)*1e3));
    subplot(3,columns,j+columns);imshow(geometricPatch(k),[0 1.2]);title('窗口内样品真值（示意）','FontName','Microsoft YaHei');
    subplot(3,columns,j+2*columns);imshow(result.images{k},[0 1.2]);title(sprintf('Measured: %d x %d',o.imageSize,o.imageSize));
end
sgtitle('上：大样品与相机物面覆盖；中：窗口对应的清晰真值；下：实际传播后采样的离焦强度', ...
    'FontName','Microsoft YaHei');
exportgraphics(f,fullfile(folder,'04_camera_acquisition.png'),'Resolution',130);
if ~o.showFigures,close(f);end

if ~o.fov.saveAcquisitionGif,return;end
f=figure('Visible','off','Position',[80 80 1350 510]);
cleanup=onCleanup(@() close(f));
for k=1:count
    subplot(1,3,1);cla;objectPanel(k);
    title('完整物面样品 + 当前相机覆盖','FontName','Microsoft YaHei');
    subplot(1,3,2);imshow(geometricPatch(k),[0 1.2]);title('当前窗口内的样品真值（示意）','FontName','Microsoft YaHei');
    subplot(1,3,3);imshow(result.images{k},[0 1.2]);title('相机记录的离焦强度','FontName','Microsoft YaHei');
    sgtitle(sprintf('Frame %d / %d | sample-to-camera distance %.2f mm',k,count,o.sampleDistances(k)*1e3));
    drawnow;rgb=print(f,'-RGBImage','-r90');[indexed,map]=rgb2ind(rgb,256);
    file=fullfile(folder,'camera_acquisition.gif');
    if k==1,imwrite(indexed,map,file,'gif','LoopCount',Inf,'DelayTime',.8);
    else,imwrite(indexed,map,file,'gif','WriteMode','append','DelayTime',.8);end
end

    function objectPanel(k)
        imagesc(coord,coord,abs(result.truth),[0 1.2]);colormap(gca,gray(256));axis image;
        hold on;
        for historyIndex=1:k
            b=edges(historyIndex,:);color=[.8 .5 .1];width=.7;
            if historyIndex==k,color=[1 0 0];width=1.8;end
            rectangle('Position',[b(1:2),b(3:4)-b(1:2)],'EdgeColor',color,'LineWidth',width);
        end
        xlabel('x [um]');ylabel('y [um]');hold off;
        xlim(viewX);ylim(viewY);
    end

    function patch=geometricPatch(k)
        % Ray projection for visualization only. Never used as measured data.
        b=result.fov.truth.planeEdges(k,:);n=o.imageSize;
        xx=b(1)+(.5:n-.5)*(b(3)-b(1))/n;
        yy=b(2)+(.5:n-.5)*(b(4)-b(2))/n;
        [xx,yy]=meshgrid(xx,yy);
        patch=interp2(abs(result.truth),xx,yy,'linear',0);
        patch(~result.cameraMasks{k})=0;
    end
end
