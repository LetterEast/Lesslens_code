function plotCoreSimulationComparison(r,folder)
%PLOTCORESIMULATIONCOMPARISON Same raw scale and evaluation region for both solvers.
idx=find(strcmp(r.names,'experimental_core'),1);if isempty(idx),return;end
visible='off';if r.options.showFigures,visible='on';end
f=figure('Visible',visible,'Name','仿真：实验核心与标定球面波对比','Position',[100 100 1350 760]);
rows=r.displayRows;cols=r.displayCols;mask=r.fov.estimated.mask;
if ~r.options.fov.displayUnion,mask=r.fov.estimated.firstMask;end
fields={r.truth,r.fields{2},r.fields{idx}};
titles={'真实样品','标定球面波：简化幅度投影','实验核心：reconstructMultiPlane'};
for k=1:3
    label=titles{k};phaseLabel='Phase [rad]';
    if k>1
        model=2;if k==3,model=idx;end
        q=r.quality(strcmp(r.quality.model,r.names{model}) & strcmp(r.quality.region,'union_fov'),:);
        label=sprintf('%s\nPSNR %.2f dB | SSIM %.3f',label,q.amplitude_psnr_db,q.amplitude_ssim);
        phaseLabel=sprintf('Phase RMSE %.4f rad | phasor SSIM %.3f',q.phase_rmse_rad,q.phase_phasor_ssim);
    end
    subplot(2,3,k);panel(abs(fields{k}),[0 1.2]);colormap(gca,gray(256));title(label,'FontName','Microsoft YaHei','Interpreter','none');
    subplot(2,3,k+3);panel(angle(fields{k}),[-1 1]);colorbar;title(phaseLabel);
end
sgtitle(sprintf('同一测量、标定与距离；%d 次迭代 | 核心 AC=%d，强度=%.3g', ...
    r.options.iterations,r.options.core.adaptiveConstraintEnabled, ...
    r.options.core.adaptiveConstraintStrength),'FontName','Microsoft YaHei');
exportgraphics(f,fullfile(folder,'06_core_vs_calibrated.png'),'Resolution',130);
if ~r.options.showFigures,close(f);end
    function panel(values,limits)
        h=imagesc(values(rows,cols),limits);h.AlphaData=double(mask(rows,cols));
        axis image;set(gca,'XTick',[],'YTick',[],'Color',[.7 .7 .7],'Box','off');
    end
end
