function [rawFcn, lb, ub, dim, fstar, info] = cec_problem(cfg, suiteName, funcId)
addpath(cfg.cec_root);
si = find(strcmp({cfg.suites.name}, suiteName), 1);
s  = cfg.suites(si);
dim = s.dim;
lb  = s.lb * ones(1, dim);
ub  = s.ub * ones(1, dim);
fstar = s.fstar(funcId);
switch suiteName
    case 'CEC2017'
        rawFcn = @(x) cec17_func(x(:), funcId);   
    case 'CEC2022'
        rawFcn = @(x) cec22_func(x(:), funcId);
    otherwise
        error('Unknown test set: %s', suiteName);
end
info = s;
end