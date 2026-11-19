function R = uav_run_one(algoKey, sc, n_ctrl, maxFEs, cfg, seed)


dim = 3*n_ctrl;
lb = zeros(1,dim); ub = cfg.space*ones(1,dim);
rawFcn = @(x) uav_cost3d(x,n_ctrl,sc.start,sc.goal,sc.obs,cfg);
MaxIter = ceil(maxFEs/cfg.N)+50;

dispName = cfg.algos{strcmp(cfg.algos(:,1),algoKey),2};
rng(seed,'twister');
t0=tic;

switch algoKey
    case {'SRDO','SRDO_NOLR'}
        bo = BudgetObjective(rawFcn,maxFEs,[],200);
        fobj = @(x) bo.eval(x);
        o = struct();
        if strcmp(algoKey,'SRDO')
        
            o.useLocalCMA=true; o.localSolver='auto'; o.localMultiK=1;
            o.localTrigger=0.8; o.localSigma=0.02; o.localFrac=0.18;
            o.localStall=999; o.localStallTrigger=0;
        end
        try
            [~,best,~] = SRDO_uni(cfg.N,MaxIter,lb,ub,dim,fobj,o,bo);
        catch ME
            if strcmp(ME.identifier,'BudgetObjective:Exhausted')
                best=bo.bestX;
            else, rethrow(ME); end
        end
        bo.finalize();
        if isempty(best), best=bo.bestX; end
        logFEs=bo.logFEs(:)'; logBest=bo.logBest(:)'; fe=bo.count;

    case {'LSHADE','IPOP_CMAES'}
        bo = BudgetObjective(rawFcn,maxFEs,[],200);
        fobj = @(x) bo.eval(x);
        try
            if strcmp(algoKey,'LSHADE')
                [~,best,~] = LSHADE_uni(cfg.N,MaxIter,lb,ub,dim,fobj,struct(),bo);
            else
                [~,best,~] = IPOP_CMAES_uni(cfg.N,MaxIter,lb,ub,dim,fobj,struct(),bo);
            end
        catch ME
            if strcmp(ME.identifier,'BudgetObjective:Exhausted')
                best=bo.bestX;
            else, rethrow(ME); end
        end
        bo.finalize();
        if isempty(best), best=bo.bestX; end
        logFEs=bo.logFEs(:)'; logBest=bo.logBest(:)'; fe=bo.count;

    otherwise
        
        fobj = rawFcn;
        switch algoKey
            case 'PSO',  [~,best,curve] = PSO(cfg.N,MaxIter,lb,ub,dim,fobj);
            case 'GA',   [~,best,curve] = GA(cfg.N,MaxIter,lb,ub,dim,fobj);
            case 'SAR',  [~,best,curve] = SAR(cfg.N,MaxIter,lb,ub,dim,fobj);
            case 'GLO',  [~,best,curve] = GLO(cfg.N,MaxIter,lb,ub,dim,fobj);
            case 'DBO',  [~,best,curve] = DBO(cfg.N,MaxIter,lb(1),ub(1),dim,fobj);
            case 'SFOA', [best,~,curve] = SFOA(cfg.N,MaxIter,lb,ub,dim,fobj);
            otherwise, error('Unknown algorithm：%s',algoKey);
        end
        curve=curve(:)'; L=numel(curve);
        logFEs=(1:L)*cfg.N; logBest=curve; fe=L*cfg.N;
end

time=toc(t0);
[~,info] = uav_cost3d(best,n_ctrl,sc.start,sc.goal,sc.obs,cfg);
R=struct('key',algoKey,'disp',dispName,'best',info.cost, ...
    'pos',best,'logFEs',logFEs,'logBest',logBest, ...
    'time',time,'fe',fe,'metrics',info);
end
