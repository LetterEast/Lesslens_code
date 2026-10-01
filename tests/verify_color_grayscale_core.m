function verify_color_grayscale_core()
% Same prepared channel, same options: verify identical saved complex fields.
setup_project;
root=fileparts(fileparts(mfilename('fullpath')));
out=fullfile(root,'outputs','color_core_verification',char(datetime('now','Format','yyyyMMdd_HHmmss_SSS')));
mkdir(out);[x,y]=meshgrid(1:64);
g=struct('M',0,'N',0,'Z',.05,'orig_M',0,'orig_N',0, ...
    'orig_size',[64 64],'padSize',[0 0],'ValidMaskHard',{{true(64),true(64)}});
inputData=struct('images',{{.5+.1*sin(x/5).*cos(y/7),.5+.08*sin(x/5).*cos(y/7)}}, ...
    'geometry',g,'wavelength',514e-9,'pixelSize',3e-6,'distanceSteps',[0 .001]);
save(fullfile(out,'channel.mat'),'inputData');
colorInput.channelFiles={'channel.mat','channel.mat','channel.mat'};
save(fullfile(out,'manifest.mat'),'colorInput');
c=loadProjectConfig('color',false);c.inputFile=fullfile(out,'manifest.mat');
c.outputRoot=fullfile(out,'color');c.reuseExistingReconstructions=true;
c.options=defaultReconstructionOptions();c.options.iterations=3;c.options.recordEvery=3;
c.options.showFigures=false;c.options.output.saveFocusPlot=false;
c.options.focus.prior=.001;c.options.focus.halfRange=0;
c.options.output.cropToValidFOV=true;c.options.output.zeroFillInvalid=true;
c.options.adaptiveConstraint.enabled=true;c.options.adaptiveConstraint.phaseMode='circular';
assert(~c.inputBrightness.enabled);
gray=struct('inputFile',fullfile(out,'channel.mat'),'options',c.options);
gray.options.output.rootFolder=fullfile(out,'gray');
[~,gr]=demo_exp_fast(gray);
cf=reconstruct_color(c);r=load(fullfile(cf,'color_fusion_result.mat'),'channelRunFolders','colorImage');
gm=jsondecode(fileread(fullfile(gr,'run_manifest.json')));gv=load(fullfile(gr,gm.latestResult),'object','field');
for k=1:3
    cr=r.channelRunFolders{k};cm=jsondecode(fileread(fullfile(cr,'run_manifest.json')));
    cv=load(fullfile(cr,cm.latestResult),'field');
    assert(max(abs(gv.field-cv.field),[],'all')<1e-12);
    channelResults=fileparts(fullfile(cr,cm.latestResult));
    grayResults=fileparts(fullfile(gr,gm.latestResult));
    assert(isequal(imread(fullfile(channelResults,'amplitude_cropped.png')), ...
        imread(fullfile(grayResults,'amplitude.png'))));
    assert(isequal(imread(fullfile(channelResults,'phase_heatmap_cropped.png')), ...
        imread(fullfile(grayResults,'phase_heatmap.png'))));
    dest=fullfile(cf,'grayscale',c.channelNames{k});
    exported=load(fullfile(dest,'reconstruction.mat'),'object');
    assert(max(abs(gv.object-exported.object),[],'all')<1e-12);
    files={'amplitude.png','phase_heatmap.png','amplitude_union.png', ...
        'phase_heatmap_union.png','originalFOV_amplitude.png','originalFOV_phase_heatmap.png'};
    for j=1:numel(files)
        assert(isequal(imread(fullfile(fileparts(fullfile(gr,gm.latestResult)),files{j})), ...
            imread(fullfile(dest,files{j}))));
    end
end
assert(all(isfinite(r.colorImage(:))));
% Fast wrapper overrides display only, including when supplied config asks
% for figures; force fresh reconstruction so this is not merely a cache test.
fast=c;fast.outputRoot=fullfile(out,'fast');fast.reuseExistingReconstructions=false;
fast.options.showFigures=true;fast.options.output.saveFocusPlot=true;
ff=reconstruct_color_fast(fast);
fr=load(fullfile(ff,'color_fusion_result.mat'),'channelRunFolders','colorImage');
assert(isequal(r.colorImage,fr.colorImage));
for k=1:3
    src=fullfile(cf,'grayscale',c.channelNames{k});dst=fullfile(ff,'grayscale',c.channelNames{k});
    for j=1:numel(files),assert(isequal(imread(fullfile(src,files{j})),imread(fullfile(dst,files{j}))));end
    saved=load(fullfile(fr.channelRunFolders{k},'run_config.mat'),'options');
    assert(~saved.options.showFigures && ~saved.options.output.saveFocusPlot);
    assert(saved.options.iterations==c.options.iterations);
end
reconstruct_color(c);r2=load(fullfile(cf,'color_fusion_result.mat'),'channelRunFolders');
assert(isequal(r.channelRunFolders,r2.channelRunFolders));
c.options.adaptiveConstraint.phaseMode='off';reconstruct_color(c);
r3=load(fullfile(cf,'color_fusion_result.mat'),'channelRunFolders');
assert(~any(strcmp(r.channelRunFolders,r3.channelRunFolders)));
if isfile(fullfile(root,'color_config.m'))
    cc=color_config();original=cc;cc.options.iterations=1;
    fresh=color_config();assert(isequaln(original,fresh)); % Independent value configuration
    assert(~contains(fileread(fullfile(root,'color_config.m')),'gray=reconstruction_config()'));
end
fprintf('Color channels equal grayscale core; cache and independent color settings passed.\n');
end
