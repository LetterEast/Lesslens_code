function field = applyAdaptiveTV(field, lambdaGradient, ...
        gradientMask, validMask, settings)
% Accelerated dual projection for spatially weighted complex TV.
dual = zeros(size(field, 1), size(field, 2), 2, 'like', field);
previousDual = dual;
lambda = cast(lambdaGradient, 'like', field);
validGradient = cast(gradientMask, 'like', field);
for subiteration = 1:settings.subiterations
    projected = dual + (1 / (8 * settings.step)) .* ...
        forwardDifference(field - settings.step .* adjointDifference(dual));
    projected = min(abs(projected), lambda) .* exp(1i .* angle(projected));
    projected = projected .* validGradient;
    momentum = subiteration / (subiteration + 3);
    dual = projected + momentum .* (projected - previousDual);
    previousDual = projected;
end
candidate = field - settings.step .* adjointDifference(previousDual);
field(validMask) = candidate(validMask);
end

function gradient = forwardDifference(field)
gradient = cat(3, field - field([2:end, end], :), ...
    field - field(:, [2:end, end]));
end

function field = adjointDifference(gradient)
vertical = gradient(:, :, 1) - gradient([end, 1:end-1], :, 1);
vertical(1, :) = gradient(1, :, 1);
vertical(end, :) = -gradient(end-1, :, 1);
horizontal = gradient(:, :, 2) - gradient(:, [end, 1:end-1], 2);
horizontal(:, 1) = gradient(:, 1, 2);
horizontal(:, end) = -gradient(:, end-1, 2);
field = vertical + horizontal;
end

