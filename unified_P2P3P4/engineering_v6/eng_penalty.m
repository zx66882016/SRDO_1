function P = eng_penalty(prob, cfg)


pMode = 1;
if strcmpi(cfg.penaltyMode,'quad'), pMode = 2; end

P.prob   = prob;
P.factor = cfg.penaltyFactor;
P.pMode  = pMode;
P.tol    = cfg.feasTol;
P.nEval  = 0;
P.histRaw  = [];
P.histViol = [];
P.histMaxg = [];

P.map   = @map_x;
P.eval  = @evaluate;
P.viol  = @violation;
P.record = @get_record;   
    function x = map_x(z)
        x = prob.lb + bsxfun(@times, z, prob.ub - prob.lb);
        mk = prob.mask;
        ds = prob.discreteSet;
        for j = 1:numel(mk)
            if ~isempty(ds) && ~isempty(ds{j})
                setv = ds{j}(:).';
                [~, k] = min(abs(x(:,j) - setv), [], 2);
                x(:,j) = setv(k);
            elseif mk(j) == 1
                x(:,j) = round(x(:,j));
            elseif mk(j) > 0
                x(:,j) = round(x(:,j)/mk(j)) * mk(j);
            end
        end
        x = min(max(x, prob.lb), prob.ub);
    end

    function [v, gnorm, gmax, graw, fraw, xraw] = violation(z)
        xraw = map_x(z);
        m = size(xraw,1);
        fraw = zeros(m,1); gnorm = zeros(m,1); gmax = zeros(m,1); graw = cell(m,1);
        for i = 1:m
            xi = xraw(i,:);
            fraw(i) = prob.rawf(xi);
            gi = prob.rawg(xi);
            graw{i} = gi(:);
            gn = max(gi(:)./prob.cScale(:), 0);
            gmax(i) = max(gn);
            if pMode == 2
                gnorm(i) = sum(gn.^2);
            else
                gnorm(i) = sum(gn);
            end
        end
        v = fraw + cfg.penaltyFactor * gnorm;
    end

    function fpen = evaluate(z)
        [fpen, vn, gm, ~, fr] = violation(z);
        if isempty(P.histRaw)
            P.histRaw  = zeros(10000,1);
            P.histViol = zeros(10000,1);
            P.histMaxg = zeros(10000,1);
        end
        m = numel(fr);
        idx = P.nEval + (1:m);
        if idx(end) > numel(P.histRaw)
            P.histRaw(2*idx(end))  = 0;
            P.histViol(2*idx(end)) = 0;
            P.histMaxg(2*idx(end)) = 0;
        end
        P.histRaw(idx)  = fr;
        P.histViol(idx) = vn;
        P.histMaxg(idx) = gm;
        P.nEval = P.nEval + m;
    end

    function [n, hr, hg, hv] = get_record()
        n  = P.nEval;
        hr = P.histRaw(1:n);
        hg = P.histMaxg(1:n);
        hv = P.histViol(1:n);
    end
end
