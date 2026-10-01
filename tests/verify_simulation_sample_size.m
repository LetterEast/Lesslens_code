function verify_simulation_sample_size()
% Changing the camera window must not resize the specimen.
setup_project;
s=struct('type','image','amplitudeFile','cameraman.tif','phaseFile','', ...
    'sizePixels',256,'offsetPixels',[0,0],'phaseScale',1,'absorptionScale',1);
a=createSimulationSample(96,96,s);b=createSimulationSample(192,96,s);
assert(isequal(a(17:272,17:272),b(65:320,65:320)));
s.sizePixels=512;failed=false;
try,createSimulationSample(96,96,s);catch e,failed=strcmp(e.identifier,'Lensless:SampleOutsideGrid');end
assert(failed);
fprintf('Specimen size stays fixed when camera window changes; oversized sample check passed.\n');
end
