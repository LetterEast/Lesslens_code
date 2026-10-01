function verify_simulation_sampling()
% Regression for the reported M=1000 carrier alias and clipped field failure.
setup_project;
o=struct('imageSize',192,'padding',96,'pixelSize',3e-6,'wavelength',514e-9, ...
    'sourceToSample',15e-3,'sourceOffsetPixels',[180,-130], ...
    'calibrationDistances',(.5:.25:1.5)*1e-3,'sampleDistances',(.5:.5:3)*1e-3);
d=checkSimulationSampling(o);assert(max(d.maxCarrierCyclesPerPixel)<.5);
o.sourceOffsetPixels=[1000,-28];assertFails(o,'Lensless:SimulationAliasing');
o.sourceOffsetPixels=[180,-130];o.calibrationDistances=[.5 5 10]*1e-3;
assertFails(o,'Lensless:CalibrationOutsideSensor');
o.calibrationDistances=[.5 1 1.5]*1e-3;o.padding=0;
assertFails(o,'Lensless:SimulationFOVClipped');
fprintf('Carrier sampling, calibration visibility and unclipped FOV checks passed.\n');
end

function assertFails(o,id)
failed=false;try,checkSimulationSampling(o);catch e,failed=strcmp(e.identifier,id);end
assert(failed,'Expected rejection: %s',id);
end
