function [fields, folder] = demo_exp_fast(config)
%DEMO_EXP_FAST Same experimental algorithm/parameters, without live figures.
root=fileparts(fileparts(mfilename('fullpath')));addpath(root);setup_project('local');
if nargin<1
    if isfile(fullfile(root,'reconstruction_config.m'))
        config=reconstruction_config();
    else
        config=loadProjectConfig('reconstruct',false);
    end
end
config.options.showFigures=false;
config.options.output.saveFocusPlot=false;
if isfield(config,'prepareInput')
    [fields,folder]=run_reconstruction(config);
else
    [fields,folder]=reconstructPreparedData(config);
end
end
