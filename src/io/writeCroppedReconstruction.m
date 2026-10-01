function writeCroppedReconstruction(object,validMask,resultsFolder)
%WRITECROPPEDRECONSTRUCTION Additional grayscale views at each saved iteration.
% Use the same valid-union crop and display scaling as standalone grayscale.
[rows,cols]=find(validMask);assert(~isempty(rows),'Empty reconstruction support.');
rows=min(rows):max(rows);cols=min(cols):max(cols);
object=object(rows,cols);mask=validMask(rows,cols);object(~mask)=0;
amplitude=normalizePercentile(abs(object),1,99,mask);
phase=phaseHeatmap(angle(object),mask);
writeMaskedPNG(amplitude,fullfile(resultsFolder,'amplitude_cropped.png'),mask);
writeMaskedPNG(phase,fullfile(resultsFolder,'phase_heatmap_cropped.png'),mask);
end
