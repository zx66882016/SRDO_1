function R = u_dispatch(cfg, suiteName, funcId, algoName, maxFEs, run, optsSrd)


if nargin < 7, optsSrd = struct(); end
here = pwd; c = onCleanup(@() cd(here));
addpath(cfg.here);
cd(cfg.cec_root); addpath(cfg.cec_root);

[rawFcn, lb, ub, dim, fstar] = cec_problem(cfg, suiteName, funcId);


suiteOff = 100000 * (find(strcmp({cfg.suites.name},suiteName),1)-1);
seed = cfg.seed_base + suiteOff + funcId*1000 + run;
rng(seed,'twister');

cps = cfg.suites(strcmp({cfg.suites.name},suiteName)).checkpoints;
cps = cps(cps <= maxFEs);
if isempty(cps), cps = maxFEs; end
if cps(end) < maxFEs, cps = [cps maxFEs]; end

bo = BudgetObjective(rawFcn, maxFEs, cps, cfg.logPoints);
fobj = @(x) bo.eval(x);
MaxIter = ceil(maxFEs/cfg.N) + 50;

t0=tic;
try
    switch algoName
        case 'SRDO'
            % Full SRDO: automatic local refinement (CMA-ES / Nelder-Mead)
            o=optsSrd; o.useLocalCMA=true; o.localSolver='auto';
            o.localMultiK=1; o.localTrigger=0.8; o.localSigma=0.02;
            o.localStall=999; o.localStallTrigger=0;
            [~,~,~] = SRDO_uni(cfg.N,MaxIter,lb,ub,dim,fobj,o,bo);
        case 'SRDO_NOLR'
            % SRDO without local refinement
            o=optsSrd;
            [~,~,~] = SRDO_uni(cfg.N,MaxIter,lb,ub,dim,fobj,o,bo);
        case 'SRDO_HP'
            o=optsSrd; o.mode='high';
            [~,~,~] = SRDO_uni(cfg.N,MaxIter,lb,ub,dim,fobj,o,bo);
        case 'AB_NM'
            % Ablation: Nelder-Mead local refinement only
            o=optsSrd; o.useLocalCMA=true; o.localSolver='nm';
            o.localTrigger=0.8; o.localStall=200; o.localStallTrigger=0;
            [~,~,~] = SRDO_uni(cfg.N,MaxIter,lb,ub,dim,fobj,o,bo);
        case 'AB_CMA'
            % Ablation: CMA-ES local refinement only
            o=optsSrd; o.useLocalCMA=true; o.localSolver='cma';
            o.localTrigger=0.8; o.localSigma=0.02; o.localStall=999; o.localStallTrigger=0;
            [~,~,~] = SRDO_uni(cfg.N,MaxIter,lb,ub,dim,fobj,o,bo);
        case 'AB_CMA_STALL'
            % Ablation: CMA-ES with stall-triggered activation
            o=optsSrd; o.useLocalCMA=true; o.localSolver='cma';
            o.localTrigger=0.8; o.localSigma=0.02; o.localStall=999;
            o.localStallTrigger=50;
            [~,~,~] = SRDO_uni(cfg.N,MaxIter,lb,ub,dim,fobj,o,bo);
        case 'AB_MULTISTART'
            % Ablation: multi-start CMA-ES (three local instances)
            o=optsSrd; o.useLocalCMA=true; o.localSolver='cma';
            o.localMultiK=3; o.localTrigger=0.8; o.localSigma=0.02;
            o.localStall=999; o.localStallTrigger=0;
            [~,~,~] = SRDO_uni(cfg.N,MaxIter,lb,ub,dim,fobj,o,bo);
        case 'AB_RANDTRIG'
            % Ablation: randomly triggered local refinement
            o=optsSrd; o.useLocalCMA=true; o.localSolver='auto';
            o.localMultiK=1; o.localTrigger=0.8; o.localSigma=0.02;
            o.localStall=999; o.localStallTrigger=0; o.localRandTrigger=true;
            [~,~,~] = SRDO_uni(cfg.N,MaxIter,lb,ub,dim,fobj,o,bo);
        case 'AB_AUTO'
            % Ablation base: automatic solver, trigger/sigma left to opts
            o=optsSrd; o.useLocalCMA=true; o.localSolver='auto';
            o.localMultiK=1; o.localStall=999; o.localStallTrigger=0;
            [~,~,~] = SRDO_uni(cfg.N,MaxIter,lb,ub,dim,fobj,o,bo);
        case 'LSHADE'
            [~,~,~] = LSHADE_uni(cfg.N,MaxIter,lb,ub,dim,fobj,struct(),bo);
        case 'IPOP_CMAES'
            [~,~,~] = IPOP_CMAES_uni(cfg.N,MaxIter,lb,ub,dim,fobj,struct(),bo);
        otherwise
            mi = floor(maxFEs/cfg.N);
            [~,~,~] = call_algorithm(algoName,cfg.N,mi,lb,ub,dim,fobj);
    end
catch ME
    if ~strcmp(ME.identifier,'BudgetObjective:Exhausted')
        rethrow(ME);
    end
end
elapsed=toc(t0);
bo.finalize();

R = struct();
R.suite=suiteName; R.func=funcId; R.algo=algoName; R.run=run; R.seed=seed;
R.maxFEs=maxFEs; R.nFEs=bo.count; R.elapsed=elapsed;
R.fstar=fstar;
R.cpFEs=bo.cpFEs;
R.cpError=max(bo.cpBest-fstar,0);
R.cpRaw=bo.cpBest;
R.logFEs=bo.logFEs;
R.logError=max(bo.logBest-fstar,0);
R.finalError=max(bo.best-fstar,0);
R.finalRaw=bo.best;
R.bestX=bo.bestX;

if isfield(optsSrd,'tag'), R.tag=optsSrd.tag; else, R.tag=algoName; end
end
