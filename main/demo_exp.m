function [fields, folder] = demo_exp(config)
%DEMO_EXP Experimental reconstruction with live progress and autofocus plots.
% Edit root reconstruction_config.m for raw images; or pass a prepared-input config.
root=fileparts(fileparts(mfilename('fullpath')));addpath(root);setup_project('local');
if nargin<1
    if isfile(fullfile(root,'reconstruction_config.m'))
        config=reconstruction_config();
    else
        config=loadProjectConfig('reconstruct',false);
    end
end
config.options.showFigures=true;
if isfield(config,'prepareInput')
    [fields,folder]=run_reconstruction(config);
else
    [fields,folder]=reconstructPreparedData(config);
end
end
