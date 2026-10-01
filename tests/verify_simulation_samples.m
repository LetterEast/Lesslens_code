function verify_simulation_samples()
% Exercise each specimen through calibration, propagation and reconstruction.
root=fileparts(fileparts(mfilename('fullpath')));addpath(root);setup_project;
out=fullfile(root,'outputs','sample_verification');if ~isfolder(out),mkdir(out);end
types={'mixed','resolution','cells','image'};
imwrite(uint8([0 255;128 64]),fullfile(out,'amplitude.png'));
imwrite(uint8([0 64;128 255]),fullfile(out,'phase.png'));
f=figure('Visible','off');cleanup=onCleanup(@() close(f));
for k=1:numel(types)
    o=struct('iterations',2,'recordEvery',1,'outputRoot',out,'sourceOffsetPixels',[30,-22]);
    o.imageSize=192;o.padding=64;o.sourceToSample=15e-3;
    o.sampleDistances=(.5:.5:3)*1e-3;
    o.fov=struct('saveAcquisitionGif',false);
    o.core=struct('enabled',false);
    o.sample=struct('type',types{k},'sizePixels',115,'offsetPixels',[0,0]);
    if strcmp(types{k},'image')
        o.sample.amplitudeFile=fullfile(out,'amplitude.png');
        o.sample.phaseFile=fullfile(out,'phase.png');
    end
    [r,folder]=demo_sim_fast(o);
    assert(all(isfinite(r.truth(:)))&&all(isfinite(r.metrics(:))));
    assert(isfile(fullfile(folder,'01_calibrated_vs_plane.png')));
    assert(strcmp(r.options.sample.type,types{k}));
    if strcmp(types{k},'cells'),assert(max(abs(abs(r.truth(:))-1))<1e-12);end
    if strcmp(types{k},'resolution'),assert(all(imag(r.truth(:))==0));end
    subplot(2,4,k,'Parent',f);imagesc(abs(r.truth(r.crop,r.crop)),[0 1]);axis image off;title(types{k});
    subplot(2,4,k+4,'Parent',f);imagesc(angle(r.truth(r.crop,r.crop)),[0 1]);axis image off;
end
exportgraphics(f,fullfile(out,'sample_gallery.png'));
% Unknown overrides must fail before expensive work.
failed=false;
try,demo_sim_fast(struct('iterationz',1));catch e,failed=contains(e.message,'Unknown simulation option');end
assert(failed);
fprintf('All four simulation sample modes and override validation passed.\n');
end
