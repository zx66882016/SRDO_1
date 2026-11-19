function [fmin, best_pos, Convergence_curve] = SAR(N, Max_Iteration, down, up, dim, fobj)
  
    social_factor = 0.8;    
    individual_factor = 0.2; 
    mu = 0.5;               

    if isscalar(down); down = down * ones(1, dim); end
    if isscalar(up); up = up * ones(1, dim); end

    X = initialization(N, dim, up, down);
    Fitness = zeros(1, N);
    for i = 1:N
        Fitness(i) = fobj(X(i, :));
    end

    [fmin, idx] = min(Fitness);
    best_pos = X(idx, :);
    global_best = best_pos;
    global_best_fit = fmin;

    Pbest = X;
    Pbest_fit = Fitness;

    Convergence_curve = inf(1, Max_Iteration);

    % Main loop
    for t = 1:Max_Iteration
        % Calculate exploration probability
        progress = t / Max_Iteration;
        exploration_prob = 1 - progress;

        for i = 1:N
            % 2.Select search strategy based on probability
            if rand < exploration_prob
               
                idx_list = randperm(N, 2);
                while any(idx_list == i)
                    idx_list = randperm(N, 2);
                end
                r1 = idx_list(1);
                r2 = idx_list(2);
             
                X_new = X(i, :) + social_factor * rand(1, dim) .* (Pbest(r1, :) - X(i, :)) ...
                                      + social_factor * rand(1, dim) .* (X(r2, :) - X(i, :));
            else
                % --- Individual search phase --
                step_size = (up - down) .* (1 - progress) * 0.1; 
                X_new = Pbest(i, :) + (rand(1, dim) - 0.5) .* step_size;
            end

            % boundary handling
            X_new = max(min(X_new, up), down);

            % Evaluate the new interpretation
            Fnew = fobj(X_new);

            % 3. Greedy selection update
            if Fnew < Fitness(i)
                X(i, :) = X_new;
                Fitness(i) = Fnew;
                
                % Update individual historical best
                if Fnew < Pbest_fit(i)
                    Pbest(i, :) = X_new;
                    Pbest_fit(i) = Fnew;
                end

                % Update global optimum
                if Fnew < global_best_fit
                    global_best = X_new;
                    global_best_fit = Fnew;
                end
            end
        end

        % 4. Abandoning the clue strategy

        for i = 1:N
            if rand < mu && Pbest_fit(i) > global_best_fit
                X(i, :) = down + (up - down) .* rand(1, dim);
                Fitness(i) = fobj(X(i, :));
                Pbest(i, :) = X(i, :);
                Pbest_fit(i) = Fitness(i);
            end
        end

        if global_best_fit < fmin
            best_pos = global_best;
            fmin = global_best_fit;
        end

        Convergence_curve(t) = global_best_fit;
    end
end

% --- helper function ---
function X = initialization(N, dim, up, down)
    X = down + (up - down) .* rand(N, dim);
end