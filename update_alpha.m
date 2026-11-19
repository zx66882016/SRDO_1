function newX = update_alpha(X, Best_X, ASF, ct, mt, lb, ub, D)

    r1 = rand();
    if r1 < 0.5
        % neighborhood search
        levy = 1 ./ (abs(randn(1, D)) .^ (1/1.5));  % Levy flight step
        a = unifrnd(-5, 5);
        b = ASF * max((1 - ct / mt), 0.4) * (ub - lb) / 30;
        newX = X + (a * b .* levy);
    else
        % approach toward the global best
        mu = unifrnd(-.1,.1, 1, D);   C1 =(10-(ct/mt))*mu;    % Eq. (7)
        newX = X + C1 .* (X - Best_X);   % Eq. (8.2
    end

end

