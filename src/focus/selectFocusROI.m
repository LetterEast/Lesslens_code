function bounds = selectFocusROI(amplitude)
%SELECTFOCUSROI Select pixel bounds [xmin ymin xmax ymax] on a focus preview.
% The caller propagates the complete field to the prior distance first.
% Closing/cancelling stops the run; it never silently falls back to background.
amplitude = double(gatherIfNeeded(amplitude));
[height,width] = size(amplitude);
fig = figure('Name','Select autofocus sample region','NumberTitle','off', ...
    'Color','w','Position',[100,80,1050,800],'Visible','on', ...
    'Tag','FocusROISelector','CloseRequestFcn',@(source,~) cancel(source));
cleanup = onCleanup(@() closeIfValid(fig));
setappdata(fig,'Confirmed',false);
setappdata(fig,'Cancelled',false);
ax = axes('Parent',fig,'Position',[0.04,0.16,0.92,0.77]);
limits = prctile(amplitude(isfinite(amplitude)),[1,99]);
if isempty(limits) || limits(2)<=limits(1), limits = []; end
imshow(amplitude,limits,'Parent',ax);
title(ax,'拖动方框到样品处，调整大小后点击“确认区域”');
side = min([256,height-1,width-1]);
roi = drawrectangle(ax,'Position',[(width-side)/2,(height-side)/2,side,side], ...
    'Color','c','DrawingArea',[0.5,0.5,width,height]);
setappdata(fig,'FocusROI',roi);
uicontrol(fig,'Style','text','Units','normalized','Position',[0.04,0.07,0.58,0.06], ...
    'BackgroundColor','w','String','仅用框内样品评价清晰度；不会裁剪最终重建图。', ...
    'FontSize',11);
uicontrol(fig,'Style','pushbutton','Units','normalized','Position',[0.66,0.06,0.16,0.06], ...
    'String','确认区域','Tag','FocusROIConfirm','Callback',@(~,~) confirm(fig));
uicontrol(fig,'Style','pushbutton','Units','normalized','Position',[0.84,0.06,0.12,0.06], ...
    'String','取消','Tag','FocusROICancel','Callback',@(~,~) cancel(fig));
if ~getappdata(fig,'Cancelled') && ~getappdata(fig,'Confirmed'), uiwait(fig); end
if ~isgraphics(fig) || ~getappdata(fig,'Confirmed') || ~isvalid(roi)
    error('Lensless:FocusSelectionCancelled','Autofocus region selection was cancelled.');
end
position = roi.Position;
bounds = [max(1,ceil(position(1))),max(1,ceil(position(2))), ...
    min(width,floor(position(1)+position(3))), ...
    min(height,floor(position(2)+position(4)))];
if bounds(3)<=bounds(1) || bounds(4)<=bounds(2)
    error('Lensless:InvalidFocusROI','Select a region at least 2 pixels wide and high.');
end
end

function confirm(fig)
setappdata(fig,'Confirmed',true);
uiresume(fig);
end

function cancel(fig)
% Keep the figure alive until setup/wait has returned, then the cleanup
% deletes it. Closing during UI construction cannot leave invalid handles.
setappdata(fig,'Confirmed',false);
setappdata(fig,'Cancelled',true);
uiresume(fig);
end

function closeIfValid(fig)
if isgraphics(fig), delete(fig); end
end
