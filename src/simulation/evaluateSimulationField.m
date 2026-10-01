function scores=evaluateSimulationField(field,truth,mask,phaseThreshold)
%EVALUATESIMULATIONFIELD Fixed-scale metrics on shared object-plane masks.
% No image normalization, clipping, registration or fitted amplitude gain.
% Phase errors are circular; low truth-amplitude pixels are excluded.
% Phase SSIM averages SSIM(cos(phi)) and SSIM(sin(phi)), each mapped to [0,1].
% It is a phasor SSIM, not ordinary grayscale SSIM on wrapped angles.
a=abs(field);ref=abs(truth);phaseMask=mask & ref>phaseThreshold;
scores.pixel_count=nnz(mask);scores.phase_pixel_count=nnz(phaseMask);
scores.amplitude_rmse=sqrt(mean((a(mask)-ref(mask)).^2));
scores.amplitude_psnr_db=20*log10(1/scores.amplitude_rmse);
scores.complex_nrmse=norm(field(mask)-truth(mask))/max(norm(truth(mask)),eps);
[~,map]=ssim(a,ref,'DynamicRange',1,'Radius',1.5);
% MATLAB uses an 11x11 Gaussian window for Radius=1.5. No masked zero fill.
valid=conv2(double(mask),ones(11),'same')==121;
scores.ssim_pixel_count=nnz(valid);
scores.amplitude_ssim=mean(map(valid));
delta=angle(field(phaseMask).*conj(truth(phaseMask)));
scores.phase_rmse_rad=sqrt(mean(delta.^2));
scores.phase_psnr_db=20*log10(2*pi/scores.phase_rmse_rad);
phaseValid=conv2(double(phaseMask),ones(11),'same')==121;
scores.phase_ssim_pixel_count=nnz(phaseValid);
p=angle(field);q=angle(truth);
[~,c]=ssim((cos(p)+1)/2,(cos(q)+1)/2,'DynamicRange',1,'Radius',1.5);
[~,s]=ssim((sin(p)+1)/2,(sin(q)+1)/2,'DynamicRange',1,'Radius',1.5);
scores.phase_phasor_ssim=mean((c(phaseValid)+s(phaseValid))/2);
% Empty regions have no score, including their normalized complex error.
if ~any(mask(:)),scores.complex_nrmse=NaN;end
end
