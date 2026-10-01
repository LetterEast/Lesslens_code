function config = loadProjectConfig(name, useLocal)
%LOADPROJECTCONFIG Load versioned defaults, then optional local overrides.
% Local scripts assign config fields and are excluded from Git. They run in
% this function's workspace; projectRoot is available for portable paths.
% A local file may override a single nested field without copying defaults.
if nargin < 2, useLocal = true; end
validNames = {'color', 'reconstruct', 'reconstruct_fast'};
name = validatestring(name, validNames);
projectRoot = fileparts(fileparts(fileparts(mfilename('fullpath'))));
config = struct();
run(fullfile(projectRoot, 'config', [name '_example.m']));
localFile = fullfile(projectRoot, 'local', 'config', [name '_local.m']);
if useLocal && isfile(localFile), run(localFile); end
end
