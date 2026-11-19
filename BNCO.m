function [fmin, best_pos, Convergence_curve] = BNCO(N, Max_Iteration, down, up, dim, fobj)
    
    if isscalar(down); down = down * ones(1, dim); end
    if isscalar(up); up = up * ones(1, dim); end


    p_explore = 0.5;      
    p_mutation = 0.1;  
    levy_beta = 1.5;    

    X = initialization(N, dim, up, down);
    Fitness = zeros(1, N);
    for i = 1:N
        Fitness(i) = fobj(X(i, :));
    end

    [fmin, idx] = min(Fitness);
    best_pos = X(idx, :);
    global_best = best_pos;
    global_best_fit = fmin;

    Convergence_curve = inf(1, Max_Iteration);

    for t = 1:Max_Iteration
        progress = t / Max_Iteration; 

        for i = 1:N
            % ---Stage 1: Foraging stage---
            if rand < (1 - progress * 0.5)  % Initially, foraging is the main activity
                % Search for prey, explore the area
                if rand < p_explore
                    % Random exploration via Lévy flight
                    step = levy_flight(dim, levy_beta);
                    scale = 0.01 + 0.1 * (1 - progress);
                    X_new = X(i, :) + scale * step .* (up - down);
                else
                    % Pouncing on prey
                    r1 = randi(N);
                    while r1 == i; r1 = randi(N); end
                    X_new = X(i, :) + rand(1, dim) .* (X(r1, :) - X(i, :)) + ...
                            0.5 * rand(1, dim) .* (global_best - X(i, :));
                end
            else
                % --- Migration Phase ---
                r2 = randi(N);
                while r2 == i; r2 = randi(N); end
                % Group migration: moving towards a weighted average of random individuals and global optimum
                X_new = X(i, :) + rand(1, dim) .* (X(r2, :) - X(i, :)) + ...
                        0.3 * rand(1, dim) .* (global_best - X(i, :));
            end

            % --- Stage 3: Survival of the fittest mechanism ---
          
            if rand < p_mutation * (1 - progress)
                % Apply Gaussian mutation to some individuals
                mutation_size = (up - down) .* (0.05 + 0.1 * rand(1, dim));
                X_new = X_new + mutation_size .* randn(1, dim);
            end

            % boundary handling
            X_new = max(min(X_new, up), down);

            % Evaluate the new interpretation
            Fnew = fobj(X_new);

            % Greedy selection
            if Fnew < Fitness(i)
                X(i, :) = X_new;
                Fitness(i) = Fnew;

                if Fnew < global_best_fit
                    global_best = X_new;
                    global_best_fit = Fnew;
                end
            end
        end

        % Update global optimum
        [current_best, idx_best] = min(Fitness);
        if current_best < global_best_fit
            global_best = X(idx_best, :);
            global_best_fit = current_best;
        end

        best_pos = global_best;
        fmin = global_best_fit;
        Convergence_curve(t) = global_best_fit;
    end
end

% ========== helper function ==========

function X = initialization(N, dim, up, down)
    X = down + (up - down) .* rand(N, dim);
end

function step = levy_flight(d, beta)

    sigma = (gamma(1 + beta) * sin(pi * beta / 2) / ...
             (gamma((1 + beta) / 2) * beta * 2^((beta - 1) / 2)))^(1 / beta);
    u = randn(1, d) * sigma;
    v = randn(1, d);
    step = u ./ abs(v).^(1 / beta);
end