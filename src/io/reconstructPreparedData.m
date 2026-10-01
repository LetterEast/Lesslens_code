function [result, runFolder] = reconstructPreparedData(config)
%RECONSTRUCTPREPAREDDATA Reconstruct a prepared MAT; return recorded fields and run folder.
% Load standard defaults and local overrides unless config is supplied.
% result 为各记录迭代的参考面复场 cell；物面图像保存在 runFolder 内。
projectRoot = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(projectRoot); setup_project;
if nargin < 1 || isempty(config), config = loadProjectConfig('reconstruct'); end
if ~isfile(config.inputFile)
    error('Lensless:MissingInput', ...
        'Prepared input MAT is missing: %s. See docs/INPUT_FORMAT.md.', config.inputFile);
end
loaded = load(config.inputFile, 'inputData');
if ~isfield(loaded, 'inputData')
    error('Lensless:InvalidInput', 'MAT file must contain inputData: %s', config.inputFile);
end
validateReconstructionInput(loaded.inputData);
[result, runFolder] = reconstructMultiPlane(loaded.inputData, config.options);
end
