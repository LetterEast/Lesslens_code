function maps = buildAdaptiveTV(coverageCount, planeUncertainty, geometry, settings)
% Adaptive object-domain TV based on overlap, disagreement and shift direction.
maps.reliableMask = coverageCount > 0;
maximumCoverage = max(coverageCount(:));
maps.confidence = coverageCount ./ max(maximumCoverage, 1);
maps.coverageRisk = (1 - maps.confidence) .^ settings.coveragePower;

uncertainty = gatherIfNeeded(real(planeUncertainty));
coreMask = maps.reliableMask & maps.confidence >= 0.8;
if any(coreMask(:))
    baseline = median(uncertainty(coreMask), 'all');
else
    baseline = median(uncertainty(maps.reliableMask), 'all');
end
uncertaintyExcess = max(uncertainty - baseline, 0);
positive = uncertaintyExcess(maps.reliableMask & uncertaintyExcess > 0);
if isempty(positive)
    uncertaintyScale = 1;
else
    uncertaintyScale = prctile(positive, 95);
end
maps.uncertainty = min(uncertaintyExcess ./ max(uncertaintyScale, eps), 1);

% Automatically protect the maximum-overlap core and increase TV outward.
distanceFromCore = bwdist(coreMask);
boundaryScale = max(distanceFromCore(maps.reliableMask), [], 'all');
maps.boundaryRisk = (distanceFromCore ./ max(boundaryScale, 1)) .^ ...
    settings.boundaryPower;
maps.boundaryRisk(~maps.reliableMask) = 0;
maps.coreProtectMask = coreMask;

maps.risk = max(settings.coverageWeight .* maps.coverageRisk, ...
    settings.uncertaintyWeight .* maps.uncertainty);
maps.risk = max(maps.risk, settings.boundaryWeight .* maps.boundaryRisk);
maps.risk(coreMask) = 0;
maps.risk(~maps.reliableMask) = 0;
maps.risk = min(max(maps.risk, 0), 1);
maps.lambdaMap = settings.lambdaMin + ...
    (settings.lambdaMax - settings.lambdaMin) .* maps.risk;
maps.lambdaMap(~maps.reliableMask) = 0;

sourceM = geometry.M; sourceN = geometry.N;
if isfield(geometry, 'orig_M'), sourceM = geometry.orig_M; end
if isfield(geometry, 'orig_N'), sourceN = geometry.orig_N; end
shiftMagnitude = max(abs([sourceM, sourceN]));
if shiftMagnitude > 0
    maps.directionWeightY = 1 + settings.directionGain * abs(sourceN) / shiftMagnitude;
    maps.directionWeightX = 1 + settings.directionGain * abs(sourceM) / shiftMagnitude;
else
    maps.directionWeightY = 1;
    maps.directionWeightX = 1;
end

verticalMask = maps.reliableMask & maps.reliableMask([2:end, end], :);
verticalMask(end, :) = false;
horizontalMask = maps.reliableMask & maps.reliableMask(:, [2:end, end]);
horizontalMask(:, end) = false;
maps.gradientMask = cat(3, verticalMask, horizontalMask);
lambdaVertical = max(maps.lambdaMap, maps.lambdaMap([2:end, end], :));
lambdaHorizontal = max(maps.lambdaMap, maps.lambdaMap(:, [2:end, end]));
lambdaVertical = settings.lambdaMin + maps.directionWeightY .* ...
    max(lambdaVertical - settings.lambdaMin, 0);
lambdaHorizontal = settings.lambdaMin + maps.directionWeightX .* ...
    max(lambdaHorizontal - settings.lambdaMin, 0);
maps.lambdaGradient = cat(3, lambdaVertical, lambdaHorizontal);
maps.lambdaGradient(~maps.gradientMask) = 0;
maps.uncertaintyBaseline = baseline;
maps.uncertaintyScale = uncertaintyScale;
maps.sourceM = sourceM;
maps.sourceN = sourceN;
end

