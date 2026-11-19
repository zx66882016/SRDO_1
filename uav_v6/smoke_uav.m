function smoke_uav()
cfg = uav_config();
cfg.runs = 2;
cfg.algos = {
    'SRDO',       'SRDO';
    'SRDO_NOLR',  'SRDO (w/o LR)';
    'LSHADE',     'L-SHADE';
    'PSO',        'PSO'};
cfg.fes24 = 3000;     
cfg.fes48 = 4000;

fprintf('\n--- A* Connectivity among four scenarios ---\n');
for s=1:numel(cfg.scenarios)
    sc=uav_scenarios(cfg.scenarios{s});
    a=uav_astar3d(sc,cfg);
    if a.found
        fprintf('%s: found, raw=%.2f (hard=%d), spline=%.2f (hard=%d), %.2fs\n',...
            sc.name,a.raw.cost,a.raw.n_hard,a.spline.cost,a.spline.n_hard,a.searchTime);
    else
        fprintf('%s: NOT FOUND\n',sc.name);
    end
end

fprintf('\n--- 24-dimensional S1 scenario 4 algorithm × 2 times ---\n');
sc=uav_scenarios('S1_sparse');
for ai=1:size(cfg.algos,1)
    for r=1:2
        seed=100000+8000+r;
        R=uav_run_one(cfg.algos{ai,1},sc,8,cfg.fes24,cfg,seed);
        fprintf('%-12s run%d cost=%10.3f feas=%d len=%.2f time=%.2fs\n',...
            cfg.algos{ai,1},r,R.best,R.metrics.feasible,R.metrics.len,R.time);
    end
end

fprintf('\n---48-dimensional S3 scene SRDO × 1 time---\n');
sc=uav_scenarios('S3_dense');
R=uav_run_one('SRDO',sc,16,cfg.fes48,cfg,300000+16000+1);
fprintf('SRDO 48 dimensions cost=%.3f feas=%d len=%.2f time=%.2fs\n',...
    R.best,R.metrics.feasible,R.metrics.len,R.time);

fprintf('\nSmoke test passed。\n');
end
