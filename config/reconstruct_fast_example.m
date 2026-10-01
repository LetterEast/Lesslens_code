% Public defaults; machine overrides belong in local/config/reconstruct_fast_local.m.
config.inputFile = fullfile(projectRoot, 'data', 'reconstruction_input.mat');
config.options = defaultReconstructionOptions('fast');
