function outputRoot = reconstruct_color(config, showFigures)
%RECONSTRUCT_COLOR Reconstruct three wavelengths; return the output directory.
% With no arguments, load color defaults and optional local overrides.
% Pass a complete config struct to run another dataset.
% 本机统一配置 color_config.m；公开调用仍支持标准通道 MAT 配置。
% 最终独立灰度图和复场保存在 outputRoot/grayscale/B、G、R。
projectRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(projectRoot); setup_project;
if nargin < 1 || isempty(config)
    if isfile(fullfile(projectRoot,'color_config.m'))
        config = color_config();
    else
        config = loadProjectConfig('color',false); % Public config is independent of legacy local scripts.
    end
end
% Fast entry shares configuration/loading and changes display flags only.
if nargin>=2
    validateattributes(showFigures,{'logical'},{'scalar'});
    config.options.showFigures=showFigures;
    if ~showFigures,config.options.output.saveFocusPlot=false;end
end
if isfield(config,'prepareInput')
    validateattributes(config.prepareInput,{'logical'},{'scalar'});
    if config.prepareInput
        setup_project('local');
        assert(exist('create_color_input_data','file')==2, ...
            'Raw color input requires the local adapter; otherwise pass prepared channel MATs.');
        acquisition = config;
        acquisition.outputFile = config.inputFile;
        config.inputFile = create_color_input_data(acquisition);
    end
end
outputRoot = reconstructColor(config);
end
