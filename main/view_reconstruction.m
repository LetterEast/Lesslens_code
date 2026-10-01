function fig = view_reconstruction(folder)
%VIEW_RECONSTRUCTION Reopen autofocus and focused images without rerunning APRW.
% Pass a run folder or a particular iteration_XXXX folder. The latest complete
% recorded iteration is selected when a run folder is supplied.
projectRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(projectRoot); setup_project;
if ~isfile(fullfile(folder, 'diagnostics', 'autofocus.mat'))
    entries = [dir(fullfile(folder, 'iteration_*', 'diagnostics', 'autofocus.mat')); dir(fullfile(folder, 'Iter_*', 'diagnostics', 'autofocus.mat'))];
    if isempty(entries)
        error('Lensless:MissingFocusResult', 'No recorded autofocus result in: %s', folder);
    end
    [~, order] = sort({entries.folder});
    folder = fileparts(entries(order(end)).folder);
end
diagnostics = load(fullfile(folder, 'diagnostics', 'autofocus.mat'));
amplitude = imread(fullfile(folder, 'results', 'amplitude.png'));
phaseRGB = imread(fullfile(folder, 'results', 'phase_heatmap.png'));
fig = plotAutofocus(diagnostics, amplitude, phaseRGB, folder);
end
