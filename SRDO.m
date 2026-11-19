function [fmin, best_pos, Convergence_curve] = SRDO(N, Max_Iteration, down, up, dim, fobj)
 
    alpha0 = 0.95;      
    p0 = 0.85;           
    elite_ratio = 0.25;  
    mutation_rate = 0.15;
    diversity_threshold = 0.02;     
    if isscalar(down); down = down * ones(1, dim); end
    if isscalar(up); up = up * ones(1, dim); end
    
    % ==================== Mixed initialization ====================
    X = stable_initialization(N, dim, up, down);
    Fitness = zeros(1, N);
    for i = 1:N
        Fitness(i) = fobj(X(i, :));
    end
    
    % ====================Elite Mechanism====================
    [best_fit, idx] = sort(Fitness);
    elite_num = max(3, round(N * elite_ratio));
    elite_pop = X(idx(1:elite_num), :);
    best_pos = elite_pop(1, :);
    fmin = best_fit(1);
    global_best = best_pos;
    global_best_fit = fmin;
    
    % ====================  Perception state initialization ====================
    Hunger = ones(N, 1);          % hunger factor
    auditory_sensitivity = 0.5 + 0.5 * rand(N, 1); % auditory acuity [0.5, 1.0]
    
    % convergence curve
    Convergence_curve = inf(1, Max_Iteration);
    last_improvement = 0;
    
    % ==================== Main loop ====================
    for t = 1:Max_Iteration
        % --- Continuous adaptive parameters
        progress = t / Max_Iteration;
        
    
        alpha = alpha0 * exp(-4 * progress^2);
        p = p0 * (1 - progress^0.8);
        mutation_factor = 0.9 * (1 - progress^0.4);
        levy_factor = 1.0 - 0.5 * progress;
        lens_k = 0.5 + 0.4 * progress;
        
        if t - last_improvement > 0.1 * Max_Iteration
            alpha = min(1.0, alpha * 1.2);
            p = min(0.95, p * 1.1);
        end
        
        % --- Update the hunger factor ---
        best_fitness = min(Fitness);
        worst_fitness = max(Fitness);
        for i = 1:N
            if Fitness(i) ~= best_fitness
                normalized_fit = (Fitness(i) - best_fitness) / (worst_fitness - best_fitness + eps);
                Hunger(i) = Hunger(i) + normalized_fit^1.5;
            end
        end
        max_hunger = max(Hunger);
        if max_hunger > 0; Hunger = Hunger / max_hunger; end
        
        % --- Unified odor guidance search & auditory localization ---
        for i = 1:N
            % 策略选择 
            rnd = rand;
            
            % ---------- Strategy A: Long-range turbulence exploration  ----------
            if rnd < 0.35  % Simulate long-distance diffusion of odor molecules
                step = stable_levy_step(dim, progress, Hunger(i), levy_factor);
                if rand < (1 - progress)  % The tendency to explore decreases as the process progresses
                    X_new = X(i, :) + alpha * step .* (up - down);
                else
                    leader = elite_pop(randi(elite_num), :);
                    X_new = X(i, :) + alpha * step .* (leader - X(i, :));
                end
                
            % ---------- Strategy B: Medium-range Perturbation Offset  ----------
            elseif rnd < 0.7  % Simulate sudden changes in wind direction or interference from companion airflows
                r1 = randi(N); r2 = randi(N);
                while r1 == i || r2 == i || r1 == r2
                    r1 = randi(N); r2 = randi(N);
                end
                F = 0.5 + 0.4 * Hunger(i) * (1 - progress); 
                X_new = X(i, :) + F * (global_best - X(i, :)) + F * (X(r1, :) - X(r2, :));
                
            % ---------- Strategy C: Short-range auditory focusing ----------
            else  % Simulate binaural effect for precise localization
                leader = elite_pop(randi(elite_num), :);
                X_new = auditory_lens_imaging(X(i, :), leader, down, up, lens_k);
            end
            
            % --- Boundary reflection processing ---
            X_new = stable_bounds(X_new, down, up);
            
            % ---Adaptive Gaussian Micro-Mutation ---
            if rand < mutation_rate * mutation_factor
                mutation_size = (up - down) .* (0.01 + 0.05 * (1 - progress) * rand(1, dim));
                X_new = X_new + mutation_size .* randn(1, dim);
                X_new = stable_bounds(X_new, down, up);
            end
            
            % ---Greedy selection and update ---
            Fnew = fobj(X_new);
            if Fnew < Fitness(i)
                X(i, :) = X_new;
                Fitness(i) = Fnew;
                if Fnew < global_best_fit
                    global_best = X_new;
                    global_best_fit = Fnew;
                    last_improvement = t;
                end
            end
        end
        
        % --- 更新精英种群 ---
        [best_fit, idx] = sort(Fitness);
        elite_pop = X(idx(1:elite_num), :);
        best_pos = elite_pop(1, :);
        fmin = best_fit(1);
        
     
        if global_best_fit < fmin
            best_pos = global_best;
            fmin = global_best_fit;
            elite_pop(1, :) = global_best;
        end
        
      
        Convergence_curve(t) = global_best_fit;
        
        % ---Dynamic maintenance of diversity---
        diversity = mean(std(X));
        if mod(t, 10) == 0 || diversity < diversity_threshold
            worst_num = max(1, round(N * 0.1));
            worst_idx = idx(end-worst_num+1:end);
            X(worst_idx, :) = stable_initialization(worst_num, dim, up, down);
            for k = 1:worst_num
                Fitness(worst_idx(k)) = fobj(X(worst_idx(k), :));
            end
            
            [best_fit, idx] = sort(Fitness);
            elite_pop = X(idx(1:elite_num), :);
        end
        
        % --- Final confirmation mechanism ---
        if t > 0.7 * Max_Iteration && mod(t, 8) == 0
            step_size = 0.02 * (up - down) * (1 - progress); 
            [refined_sol, refined_fit] = pattern_search(global_best, down, up, fobj, step_size, 8);
            if refined_fit < global_best_fit
                global_best = refined_sol;
                global_best_fit = refined_fit;
                elite_pop(1, :) = refined_sol;
                Fitness(idx(1)) = refined_fit;
                Convergence_curve(t) = global_best_fit; 
            end
        end
    end
    
 
    [best_pos, fmin] = pattern_search(global_best, down, up, fobj, 0.001*(up-down), 12);
    if fmin < global_best_fit
        global_best_fit = fmin;
        best_pos = best_pos;
    else
        fmin = global_best_fit;
        best_pos = global_best;
    end
    Convergence_curve(end) = fmin;
end

% ========================================================================
%                    helper function
% ========================================================================

function X = stable_initialization(N, dim, up, down)
    
    X = zeros(N, dim);
    chaos_num = round(N * 0.6);
    for i = 1:chaos_num
        chaos = rand(1, dim);
        for k = 1:8 
            chaos = tent_map(chaos);
        end
        X(i, :) = down + (up - down) .* chaos;
    end
    rand_num = N - chaos_num;
    if rand_num > 0
        X(chaos_num+1:end, :) = down + (up - down) .* rand(rand_num, dim);
    end
end

function seq = tent_map(seq)
    
    a = 0.7;
    for j = 1:length(seq)
        if seq(j) < a
            seq(j) = seq(j) / a;
        else
            seq(j) = (1 - seq(j)) / (1 - a);
        end
    end
end

function vec = stable_levy_step(d, progress, hunger_factor, levy_factor)
    
    beta = 1.9 - 0.7 * progress + 0.1 * hunger_factor;
    beta = max(1.2, min(beta, 1.95));
    
    sigma = (gamma(1 + beta) * sin(pi * beta / 2) / ...
            (gamma((1 + beta) / 2) * beta * 2^((beta - 1) / 2)))^(1 / beta);
    u = randn(1, d) * sigma;
    v = randn(1, d);
    v(v == 0) = 1e-10; 
    step = u ./ abs(v).^(1 / beta);
    

    scale = 0.01 + 0.2 * (1 - progress)^1.5;
    vec = scale * step * levy_factor;
end

function X_new = auditory_lens_imaging(X, elite, lb, ub, k)
    
    dynamic_lb = min(X, elite);
    dynamic_ub = max(X, elite);
    
    
    range = dynamic_ub - dynamic_lb;
    if any(range < 1e-12)
        dynamic_lb = max(dynamic_lb - 0.5, lb);
        dynamic_ub = min(dynamic_ub + 0.5, ub);
    end
    
    dynamic_lb = max(dynamic_lb, lb);
    dynamic_ub = min(dynamic_ub, ub);
    
 
    center = (dynamic_lb + dynamic_ub) / 2;
    lens_point = center + (center - X) / k;
    
   
    X_new = 0.6 * lens_point + 0.4 * elite;
end

function X = stable_bounds(X, Lb, Ub)
  
    below = X < Lb;
    X(below) = 2 * Lb(below) - X(below);
    
    above = X > Ub;
    X(above) = 2 * Ub(above) - X(above);
    
    
    X = max(min(X, Ub), Lb);
    
  
    X(X < Lb) = Lb(X < Lb);
    X(X > Ub) = Ub(X > Ub);
end

function [pattern_sol, pattern_fit] = pattern_search(sol, lb, ub, fobj, step_size, max_iter)
    
    if nargin < 6
        max_iter = 12;
    end
    n = length(sol);
    current_sol = sol;
    current_fit = fobj(sol);
    
    for iter = 1:max_iter
        improved = false;
        for i = 1:n
            candidate = current_sol;
            candidate(i) = min(current_sol(i) + step_size(i), ub(i));
            candidate_fit = fobj(candidate);
            if candidate_fit < current_fit
                current_sol = candidate;
                current_fit = candidate_fit;
                improved = true;
                continue;
            end
            % 反向探测
            candidate(i) = max(current_sol(i) - step_size(i), lb(i));
            candidate_fit = fobj(candidate);
            if candidate_fit < current_fit
                current_sol = candidate;
                current_fit = candidate_fit;
                improved = true;
            end
        end
        if ~improved
            step_size = 0.5 * step_size;
            if max(step_size) < 1e-10
                break;
            end
        end
    end
    pattern_sol = current_sol;
    pattern_fit = current_fit;
end