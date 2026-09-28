function [weights, coverageCount] = coverageFusionWeights(hardMasks)
%COVERAGEFUSIONWEIGHTS Equal participation per locally covered pixel.
% Sensor edges participate fully. No support means zero measurement weight.
support = double(cat(3, hardMasks{:}) ~= 0);
coverageCount = sum(support, 3);
weights = support ./ max(coverageCount, 1);
end
