function verify_manual_focus_roi()
%VERIFY_MANUAL_FOCUS_ROI Automated graphics test; requires figure support.
% Explicitly run this test to exercise confirmation/cancellation and reuse.
projectRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(projectRoot);
setup_project;
selectionCount = 0;
cancelSelection = false;
driver = timer('ExecutionMode','fixedSpacing','Period',0.25, ...
    'TasksToExecute',160,'TimerFcn',@driveSelection,'StopFcn',@closeSelector);
cleanup = onCleanup(@() disposeTimer(driver));
start(driver);
input = syntheticInput();
options = defaultReconstructionOptions('standard');
options.iterations = 2;
options.recordEvery = 1;
options.focus.method = 'adfrft';
options.focus.roiMode = 'manual';
options.focus.prior = 1e-3;
options.focus.halfRange = 0;
options.tv.enabled = false;
options.showFigures = false;
options.output.rootFolder = fullfile(projectRoot,'outputs','manual_roi_verification', ...
    char(datetime('now','Format','yyyyMMdd_HHmmss_SSS')));
[~,folder] = reconstructMultiPlane(input,options);
assert(selectionCount==1,'Two recorded iterations must require only one selection.');
for iteration = 1:2
    focus = load(fullfile(folder,sprintf('iteration_%04d',iteration),'diagnostics','autofocus.mat'));
    assert(isequal(focus.roiBounds,[10,12,25,22]));
    assert(strcmp(focus.roiMode,'manual'));
end
fig = view_reconstruction(folder);
assert(~isempty(findobj(fig,'Tag','FocusROIOutline')));
close(fig);
cancelSelection = true;
cancelled = false;
try
    selectFocusROI(ones(24,32));
catch exception
    cancelled = strcmp(exception.identifier,'Lensless:FocusSelectionCancelled');
end
assert(cancelled);
assert(isempty(findall(groot,'Tag','FocusROISelector')));
clear cleanup;
fprintf('Manual ROI confirmation, reuse, saved overlay and cancellation passed.\n');

    function driveSelection(~,~)
        selector = findall(groot,'Type','figure','Tag','FocusROISelector');
        if isempty(selector) || isempty(getappdata(selector(1),'FocusROI')), return; end
        selector = selector(1);
        if ~strcmp(selector.WaitStatus,'waiting'), return; end
        if isequal(getappdata(selector,'TestHandled'),true), return; end
        setappdata(selector,'TestHandled',true);
        if cancelSelection
            close(selector); % Exercise the actual window-close callback.
            return;
        end
        roi = getappdata(selector,'FocusROI');
        roi.Position = [10,12,15,10];
        selectionCount = selectionCount+1;
        button = findobj(selector,'Tag','FocusROIConfirm');
        callback = button.Callback;
        callback(button,[]);
    end

    function closeSelector(~,~)
        delete(findall(groot,'Type','figure','Tag','FocusROISelector'));
    end
end

function disposeTimer(driver)
if isvalid(driver), stop(driver); delete(driver); end
end
