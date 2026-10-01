function verify_simulation_metrics()
% Metrics must preserve radiometric error and handle phase wrapping / empty FOV.
setup_project;
truth=ones(40);mask=false(40);mask(8:32,8:32)=true;
q=evaluateSimulationField(truth,truth,mask,.05);
assert(q.amplitude_rmse==0 && isinf(q.amplitude_psnr_db));
assert(abs(q.amplitude_ssim-1)<1e-12 && abs(q.phase_phasor_ssim-1)<1e-12);
q=evaluateSimulationField(.5*truth,truth,mask,.05);
assert(abs(q.amplitude_psnr_db-20*log10(2))<1e-10);
q=evaluateSimulationField(exp(1i*(-pi+.01))*truth,exp(1i*(pi-.01))*truth,mask,.05);
assert(abs(q.phase_rmse_rad-.02)<1e-10);
q=evaluateSimulationField(truth,truth,false(40),.05);
assert(q.pixel_count==0 && isnan(q.amplitude_psnr_db) && isnan(q.amplitude_ssim));
q=evaluateSimulationField(truth,zeros(40),mask,.05);
assert(q.phase_pixel_count==0 && isnan(q.phase_rmse_rad));
% SSIM windows may not leak outside a shared region.
changed=truth;changed(~mask)=100;
q=evaluateSimulationField(changed,truth,mask,.05);
assert(abs(q.amplitude_ssim-1)<1e-9);
fprintf('Metric scale, circular phase, empty regions and SSIM window masking passed.\n');
end
