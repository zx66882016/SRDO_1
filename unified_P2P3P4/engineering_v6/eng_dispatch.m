function R = eng_dispatch(cfg, pid, algoName, run)

prob = cfg.problems(pid);
P = eng_penalty(prob, cfg);
D = prob.dim;

seed = cfg.seed_base + cfg.eng_offset + pid*1000 + run;
rng(seed,'twister');

bo = BudgetObjective(@(z) P.eval(z), cfg.maxFEs, cfg.maxFEs, cfg.logPoints);
fobj = @(x) bo.eval(x);
lbN = zeros(1,D); ubN = ones(1,D);
MaxIter = ceil(cfg.maxFEs/cfg.N) + 50;

t0 = tic;
try
    switch algoName
        case 'SRDO'
            % Full SRDO: automatic local refinement
            o=struct(); o.useLocalCMA=true; o.localSolver='auto';
            o.localMultiK=1; o.localTrigger=0.8; o.localSigma=0.02;
            o.localStall=999; o.localStallTrigger=0;
            SRDO_uni(cfg.N,MaxIter,lbN,ubN,D,fobj,o,bo);
        case 'SRDO_NOLR'
            % SRDO without local refinement
            SRDO_uni(cfg.N,MaxIter,lbN,ubN,D,fobj,struct(),bo);
        case 'SRDO_HP'
            o=struct(); o.mode='high';
            SRDO_uni(cfg.N,MaxIter,lbN,ubN,D,fobj,o,bo);
        case 'LSHADE'
            LSHADE_uni(cfg.N,MaxIter,lbN,ubN,D,fobj,struct(),bo);
        case 'IPOP_CMAES'
            IPOP_CMAES_uni(cfg.N,MaxIter,lbN,ubN,D,fobj,struct(),bo);
        otherwise
            mi = floor(cfg.maxFEs/cfg.N);
            call_algorithm(algoName,cfg.N,mi,lbN,ubN,D,fobj);
    end
catch ME
    if ~strcmp(ME.identifier,'BudgetObjective:Exhausted')
        rethrow(ME);
    end
end
elapsed = toc(t0);
bo.finalize();

[n, histRaw, histMaxg] = P.record();
feasFlag = histMaxg <= cfg.feasTol;

bestFeas = nan(n,1);
runBest = inf;
for i = 1:n
    if feasFlag(i) && histRaw(i) < runBest
        runBest = histRaw(i);
    end
    if isfinite(runBest), bestFeas(i) = runBest; end
end
logBestFeas = bestFeas(bo.logFEs);

xraw = P.map(bo.bestX);
fraw = prob.rawf(xraw);
gval = prob.rawg(xraw);
gnorm = max(gval(:)./prob.cScale(:),0);
maxViol = max(gnorm);
if strcmpi(cfg.penaltyMode,'quad')
    viol = sum(gnorm.^2);
else
    viol = sum(gnorm);
end
feasible = maxViol <= cfg.feasTol;

R = struct();
R.pid=pid; R.tag=prob.tag; R.algo=algoName; R.run=run; R.seed=seed;
R.maxFEs=cfg.maxFEs; R.nFEs=bo.count; R.elapsed=elapsed;
R.feasible=feasible; R.maxViol=maxViol; R.viol=viol;
R.finalRaw = fraw + 0;
if ~feasible, R.finalRaw = NaN; end   
R.finalRawAll = fraw;                 
R.bestXraw = xraw; R.g = gval(:).';
R.logFEs = bo.logFEs; R.logBestFeas = logBestFeas;
end
