function verify_constraint_phase_wrap()
% Nearly identical complex fields across pi must not acquire a finite seam.
setup_project;
n=64;[x,y]=meshgrid((1:n)-floor(n/2)-1);disk=x.^2+y.^2<12^2;
g=struct('Z',.05,'M',0,'N',0,'orig_M',0,'orig_N',0, ...
    'padSize',[0 0],'orig_size',[n n],'ValidMaskHard',{{true(n)}});
d=struct('geometry',g,'images',{{ones(n)}},'distanceSteps',0, ...
    'pixelSize',3e-6,'wavelength',514e-9);
o=defaultReconstructionOptions();s=o.adaptiveConstraint;
s.distance=0;s.whiteAmplitude=1;s.edgeWidth=4;s.strength=.05;
illum=exp(1i*2*pi/d.wavelength*sqrt((x*d.pixelSize).^2+(y*d.pixelSize).^2+g.Z^2));
a=ones(n);a(disk)=.2;delta=1e-6;
u1=a.*exp(1i*(pi-delta)*disk).*illum;
u2=a.*exp(1i*(-pi+delta)*disk).*illum;
modes={'wrapped','circular','off'};difference=zeros(1,3);
for k=1:3
    s.phaseMode=modes{k};v1=applyAdaptiveConstraint(u1,d,s);
    v2=applyAdaptiveConstraint(u2,d,s);
    difference(k)=max(abs(v1(disk)-v2(disk)));
    assert(all(isfinite(v1(:))) && all(isfinite(v2(:))));
    if k==3
        assert(max(abs(angle(v1(:).*conj(u1(:)))))<1e-12);
    end
end
assert(difference(1)>1e-3);
assert(difference(2)<1e-5 && difference(3)<1e-5);
s.strength=0;assert(isequal(applyAdaptiveConstraint(u1,d,s),u1));
fprintf('Across pi: wrapped %.6g, circular %.6g, off %.6g complex-field difference.\n',difference);
end
