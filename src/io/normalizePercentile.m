function image = normalizePercentile(data, lowPercentile, highPercentile, mask)
%NORMALIZEPERCENTILE Shared grayscale display scaling for every channel.
data = double(gatherIfNeeded(data));
if nargin<4, mask = true(size(data)); end
limits = prctile(data(mask), [lowPercentile, highPercentile]);
image = (data - limits(1)) / max(limits(2) - limits(1), eps);
image = min(max(image, 0), 1);
end
