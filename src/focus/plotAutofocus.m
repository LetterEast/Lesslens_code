function fig = plotAutofocus(diagnostics, amplitude, phaseRGB, label, fig, visible)
%PLOTAUTOFOCUS Plot the sampled metric and the saved sample-plane outputs.
% Reuse one figure across recorded iterations. Supports older autofocus MAT
% files that only contain distances, focusMetric and bestDistance.
if nargin < 5, fig = []; end
if nargin < 6, visible = 'on'; end
if isempty(fig) || ~isgraphics(fig)
    fig = figure('Name', 'APRW autofocus and focused result', ...
        'NumberTitle', 'off', 'Color', 'w', 'Visible', visible, ...
        'Position', [100, 100, 1320, 480]);
else
    clf(fig);
    set(fig, 'Visible', visible);
end
layout = tiledlayout(fig, 1, 3, 'Padding', 'compact', 'TileSpacing', 'compact');
ax = nexttile(layout);
distances = double(diagnostics.distances(:));
metric = double(diagnostics.focusMetric(:));
[~, index] = min(abs(distances-diagnostics.bestDistance));
plot(ax, distances*1e3, metric, '-o', 'MarkerSize', 3, 'LineWidth', 1.3, ...
    'Tag', 'FocusCurve');
hold(ax, 'on');
plot(ax, distances(index)*1e3, metric(index), 'rp', 'MarkerSize', 12, ...
    'MarkerFaceColor', 'r', 'Tag', 'FocusBest');
xline(ax, diagnostics.bestDistance*1e3, '--r', 'HandleVisibility', 'off');
if isfield(diagnostics, 'prior')
    xline(ax, diagnostics.prior*1e3, ':k', 'Prior', 'HandleVisibility', 'off');
end
xlabel(ax, 'Sample distance (mm)');
ylabel(ax, 'Focus metric (larger is sharper)');
grid(ax, 'on');
span = max(metric)-min(metric);
margin = max(span*0.08, max(1, max(abs(metric)))*1e-6);
ylim(ax, [min(metric)-margin, max(metric)+margin]);
status = 'Maximum sampled sharpness';
if numel(distances) == 1
    status = 'Single distance; no focus scan';
elseif max(metric)-min(metric) <= eps(max(1, max(abs(metric))))*16
    status = 'Flat curve; focus is not resolved';
elseif index == 1 || index == numel(distances)
    status = 'Peak at search boundary; review scan range';
end
methodLabel = 'FFT high-pass amplitude';
if isfield(diagnostics,'methodLabel'), methodLabel = diagnostics.methodLabel; end
title(ax, {methodLabel, sprintf('Selected: %.4f mm', diagnostics.bestDistance*1e3), status});
legend(ax, {'Sampled metric', 'Selected distance'}, 'Location', 'best');
ax = nexttile(layout);
imshow(amplitude, [], 'Parent', ax);
title(ax, 'Focused amplitude (saved output)');
if isfield(diagnostics,'roiMode') && strcmp(diagnostics.roiMode,'manual') && ...
        isfield(diagnostics,'displayBounds')
    b = diagnostics.displayBounds; r = diagnostics.roiBounds;
    r = [max(r(1),b(1)),max(r(2),b(2)),min(r(3),b(3)),min(r(4),b(4))];
    if r(3)>r(1) && r(4)>r(2)
        rectangle(ax,'Position',[r(1)-b(1)+1,r(2)-b(2)+1,r(3)-r(1),r(4)-r(2)], ...
            'EdgeColor','c','LineWidth',1.2,'Tag','FocusROIOutline');
    end
    title(ax,'Focused amplitude | cyan: manual focus ROI');
end
ax = nexttile(layout);
imshow(phaseRGB, 'Parent', ax);
title(ax, 'Focused phase (saved output)');
title(layout, label, 'Interpreter', 'none');
drawnow;
end
