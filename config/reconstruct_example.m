% Public defaults; machine overrides belong in local/config/reconstruct_local.m.
config.inputFile = fullfile(projectRoot, 'data', 'reconstruction_input.mat');
config.options = defaultReconstructionOptions('standard');
% Optional experimental constraint; see docs/RECONSTRUCTION.md.
% config.options.adaptiveConstraint.enabled = true;
% config.options.adaptiveConstraint.strength = 0.05;
% config.options.adaptiveConstraint.distance = []; % autofocus once, or known z [m]
