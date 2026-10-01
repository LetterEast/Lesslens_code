function exportGrayscaleResult(r,folder,output,sourceResult)
%EXPORTGRAYSCALERESULT Render an existing full channel field like grayscale.
% No propagation, autofocus or TV is repeated. This is before RGB alignment,
% gain/colour correction, common-field cropping and fusion.
if ~isfolder(folder),mkdir(folder);end
validMask=logical(r.validMask);object=r.object;
crop=true;if isfield(output,'cropToValidFOV'),crop=output.cropToValidFOV;end
zero=true;if isfield(output,'zeroFillInvalid'),zero=output.zeroFillInvalid;end
b=[1 1 size(object,2) size(object,1)];
if crop
    [rows,cols]=find(validMask);assert(~isempty(rows),'Empty grayscale field of view.');
    b=[min(cols) min(rows) max(cols) max(rows)];
end
rows=b(2):b(4);cols=b(1):b(3);
object=object(rows,cols);validMask=validMask(rows,cols);
objectOriginalFOV=r.objectOriginalFOV;originalValidMask=r.originalValidMask;
if zero
    object(~validMask)=0;objectOriginalFOV(~originalValidMask)=0;
end
amplitudeUnion=normalizePercentile(abs(object),1,99,validMask);
[phaseUnion,phaseDisplayLimits]=phaseHeatmap(angle(object),validMask);
amplitudeOriginal=normalizePercentile(abs(objectOriginalFOV),1,99);
[phaseOriginal,originalPhaseDisplayLimits]=phaseHeatmap(angle(objectOriginalFOV));
first=isfield(output,'defaultToOriginalFOV') && output.defaultToOriginalFOV;
if first
    amplitude=amplitudeOriginal;phase=phaseOriginal;displayMask=originalValidMask;
    defaultView='first_camera_sample_fov';
else
    amplitude=amplitudeUnion;phase=phaseUnion;displayMask=validMask;
    defaultView='geometric_sample_union';
end
writeMaskedPNG(amplitude,fullfile(folder,'amplitude.png'),displayMask);
writeMaskedPNG(phase,fullfile(folder,'phase_heatmap.png'),displayMask);
writeMaskedPNG(amplitudeUnion,fullfile(folder,'amplitude_union.png'),validMask);
writeMaskedPNG(phaseUnion,fullfile(folder,'phase_heatmap_union.png'),validMask);
imwrite(amplitudeOriginal,fullfile(folder,'originalFOV_amplitude.png'));
imwrite(phaseOriginal,fullfile(folder,'originalFOV_phase_heatmap.png'));
% Preserve numerical fields and provenance; avoid duplicating full detector field.
bounds=b+[r.bounds(1)-1 r.bounds(2)-1 r.bounds(1)-1 r.bounds(2)-1];
originalBounds=r.originalBounds;focusDistance=r.focusDistance;
sampleCoverage=r.sampleCoverage(rows,cols);trustedMask=r.trustedMask(rows,cols);
propagationMargin=r.propagationMargin;
save(fullfile(folder,'reconstruction.mat'),'object','validMask','objectOriginalFOV', ...
    'originalValidMask','bounds','originalBounds','focusDistance','sampleCoverage', ...
    'trustedMask','propagationMargin','phaseDisplayLimits','originalPhaseDisplayLimits', ...
    'defaultView','sourceResult','output');
end
