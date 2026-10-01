function [rgb, limits] = phaseHeatmap(phase, mask)
%PHASEHEATMAP Shared phase display; robust limits from the central field.
phase = double(gatherIfNeeded(phase));
if nargin<2, mask = true(size(phase)); end
rowMargin = floor(size(phase, 1) * 0.25);
columnMargin = floor(size(phase, 2) * 0.25);
central = phase(rowMargin+1:end-rowMargin,columnMargin+1:end-columnMargin);
centralMask = mask(rowMargin+1:end-rowMargin,columnMargin+1:end-columnMargin);
values = central(centralMask);
if isempty(values), values = phase(mask); end
limits = prctile(values, [1, 99]);
if limits(2) <= limits(1), limits(2) = limits(1) + eps; end
normalized = (phase - limits(1)) / (limits(2) - limits(1));
normalized = min(max(normalized, 0), 1);
indices = round(normalized * 255) + 1;
rgb = ind2rgb(indices, hot(256));
end
