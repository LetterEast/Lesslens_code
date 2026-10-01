function [result, runFolder] = reconstruct_fast(config)
%RECONSTRUCT_FAST Use the headless profile through the shared gray entry point.
projectRoot = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(projectRoot); setup_project;
if nargin < 1 || isempty(config), config = loadProjectConfig('reconstruct_fast'); end
[result, runFolder] = reconstruct(config);
end
