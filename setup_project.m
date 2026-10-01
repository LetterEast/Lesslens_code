function setup_project(mode)
%SETUP_PROJECT Add public code; use setup_project('local') for private tools.
% Tests, datasets and results are never recursively added to the path.
if nargin < 1, mode = 'core'; end
mode = validatestring(mode, {'core','local'});
projectRoot = fileparts(mfilename('fullpath'));
addpath(fullfile(projectRoot, 'main'), genpath(fullfile(projectRoot, 'src')), ...
    fullfile(projectRoot, 'examples'));
if strcmp(mode, 'local')
    folders = {fullfile(projectRoot,'local','src'), fullfile(projectRoot,'local','tools'), fullfile(projectRoot,'local','experiments')};
    for index = 1:numel(folders)
        if isfolder(folders{index}), addpath(folders{index}); end
    end
end
end
