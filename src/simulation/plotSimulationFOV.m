function plotSimulationFOV(result,folder)
%PLOTSIMULATIONFOV Compare first-frame crop and union of the SAME reconstruction.
% Geometric footprints indicate coverage, not a guarantee of recovered detail.
fov=result.fov;support=fov.truth;t=result.fields{2};truth=result.truth;
valid=support.mask|fov.estimated.mask;[rr,cc]=find(valid);
rows=max(1,min(rr)-5):min(size(t,1),max(rr)+5);
cols=max(1,min(cc)-5):min(size(t,2),max(cc)+5);
visible='off';if result.options.showFigures,visible='on';end
f=figure('Visible',visible,'Name','仿真：首帧与联合视场','Position',[100 100 1350 760]);
subplot(2,3,1);panel(abs(truth),support.mask,[0 1.2]);colormap(gca,gray(256));title('真值：多帧联合视场');
subplot(2,3,2);panel(abs(t),fov.estimated.firstMask,[0 1.2]);colormap(gca,gray(256));title('同一多帧重建：仅首帧视场');
subplot(2,3,3);panel(abs(t),fov.estimated.mask,[0 1.2]);colormap(gca,gray(256));title('标定球面波重建：联合视场');
subplot(2,3,4);panel(double(support.coverage),support.mask,[0 numel(result.images)]);colorbar;title('真值几何覆盖帧数');
subplot(2,3,5);panel(abs(abs(t)-abs(truth)),fov.addedMask,[0 .5]);colorbar;title('新增区域：幅度绝对误差');
phaseMask=fov.addedMask & abs(truth)>result.options.evaluation.phaseAmplitudeThreshold;
subplot(2,3,6);panel(abs(angle(t.*conj(truth))),phaseMask,[0 pi]);colorbar;title('新增区域：圆周相位误差 [rad]');
sgtitle(sprintf('几何覆盖增加 %.2f%%；样品范围内有效面积增加 %.2f%%\n青线：首帧，黄线：联合视场；超出样品的部分为零', ...
    100*(fov.areaGain-1),100*(fov.sampleAreaGain-1)),'FontName','Microsoft YaHei');
exportgraphics(f,fullfile(folder,'02_expanded_fov.png'),'Resolution',130);
if ~result.options.showFigures,close(f);end

q=result.quality;regions={'first_fov','union_fov','added_fov','overlap_fov'};
metrics={'amplitude_psnr_db','amplitude_ssim','phase_rmse_rad','phase_phasor_ssim'};
titles={'Amplitude PSNR [dB] (higher is better)','Amplitude SSIM (higher is better)', ...
    'Circular phase RMSE [rad] (lower is better)','Phase phasor SSIM (higher is better)'};
f=figure('Visible',visible,'Name','仿真：客观评价','Position',[120 120 1150 680]);
for k=1:4
    values=nan(4,numel(result.names));
    for j=1:4
        for m=1:numel(result.names)
            row=strcmp(q.model,result.names{m}) & strcmp(q.region,regions{j});
            values(j,m)=q.(metrics{k})(row);
        end
    end
    subplot(2,2,k);bar(values);xticks(1:4);xticklabels({'First','Union','Added','Overlap'});
    title(titles{k});grid on;
    if k==1,legend(result.names,'Interpreter','none','Location','best');end
end
sgtitle('Shared truth masks; no display normalization; empty/narrow regions may have NaN SSIM');
exportgraphics(f,fullfile(folder,'03_quality_metrics.png'),'Resolution',130);
if ~result.options.showFigures,close(f);end

    function panel(values,mask,limits)
        % Physical coordinates in the object plane, identical for every panel.
        px=result.options.pixelSize*1e6;
        xx=(cols-floor(size(t,2)/2)-1)*px;yy=(rows-floor(size(t,1)/2)-1)*px;
        im=imagesc(xx,yy,values(rows,cols),limits);im.AlphaData=double(mask(rows,cols));
        set(gca,'Color',[.65 .65 .65],'FontName','Microsoft YaHei');axis image;
        hold on;contour(xx,yy,double(support.firstMask(rows,cols)),[.5 .5],'c','LineWidth',.8);
        contour(xx,yy,double(support.mask(rows,cols)),[.5 .5],'y','LineWidth',.8);
        xlabel('x [um]');ylabel('y [um]');
    end
end
