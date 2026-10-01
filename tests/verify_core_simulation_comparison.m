function verify_core_simulation_comparison()
setup_project;
o=struct('imageSize',192,'padding',64,'sourceToSample',15e-3, ...
    'sourceOffsetPixels',[30,-22],'sampleDistances',(.5:.5:3)*1e-3, ...
    'iterations',6,'recordEvery',6);
o.sample=struct('type','mixed');o.fov=struct('saveAcquisitionGif',false);
o.core=struct('enabled',true);
[r,folder]=demo_sim_fast(o);
assert(numel(r.fields)==4 && height(r.quality)==16);
assert(strcmp(r.names{4},'experimental_core'));
assert(r.coreInfo.iterationsCompleted==o.iterations);
assert(r.coreInfo.settings.tv.enabled && r.coreInfo.settings.adaptiveConstraint.enabled);
assert(r.coreInfo.settings.focus.halfRange==0);
assert(all(isfinite(r.fields{4}(:))));
assert(isfile(fullfile(folder,'06_core_vs_calibrated.png')));
assert(isfile(fullfile(r.coreInfo.runFolder,'run_manifest.json')));
assert(isequal(size(r.histories),[o.iterations 4]));
fprintf('Real experimental core invoked; four-arm comparison and common scoring passed.\n');
end
