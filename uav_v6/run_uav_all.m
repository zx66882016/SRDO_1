function run_uav_all()


cfg = uav_config();
nA = size(cfg.algos,1);
nS = numel(cfg.scenarios);


outFile = fullfile(cfg.data_dir,'uav_results.mat');
if exist(outFile,'file')
    load(outFile,'data24','data48','astar','done24','done48','doneAstar');
else
    data24 = cell(nS,nA,cfg.runs);
    data48 = cell(nA,cfg.runs);
    astar  = cell(nS,1);
    done24=false(nS,nA,cfg.runs);
    done48=false(nA,cfg.runs);
    doneAstar=false(nS,1);
end

for s = 1:nS
    sc = uav_scenarios(cfg.scenarios{s});
    for ai = 1:nA
        key = cfg.algos{ai,1};
        for r = 1:cfg.runs
            if done24(s,ai,r), continue; end
            seed = s*100000 + 8*1000 + r;  
            data24{s,ai,r} = uav_run_one(key,sc,8,cfg.fes24,cfg,seed);
            done24(s,ai,r)=true;
        end
        fprintf('[24 dimensions] scene%d/4 %s Complete\n',s,key);
        save(outFile,'data24','data48','astar','done24','done48','doneAstar');
    end
end

%% (2)48-dimensional high-dimensional extension (only for S3 dense scenarios)
s3 = find(strcmp(cfg.scenarios,cfg.hd_scenario));
sc = uav_scenarios(cfg.hd_scenario);
for ai = 1:nA
    key = cfg.algos{ai,1};
    for r = 1:cfg.runs
        if done48(ai,r), continue; end
        seed = s3*100000 + 16*1000 + r;
        data48{ai,r} = uav_run_one(key,sc,16,cfg.fes48,cfg,seed);
        done48(ai,r)=true;
    end
    fprintf('[48-dimensional] %s Complete\n',key);
    save(outFile,'data24','data48','astar','done24','done48','doneAstar');
end

%% (3) A* Deterministic Baseline
for s = 1:nS
    if doneAstar(s), continue; end
    sc = uav_scenarios(cfg.scenarios{s});
    astar{s} = uav_astar3d(sc,cfg);
    doneAstar(s)=true;
    if astar{s}.found
        fprintf('[A*] Scene %d found a path, search  %.2f s, raw cost %.1f, spline cost %.1f\n', ...
            s,astar{s}.searchTime,astar{s}.raw.cost,astar{s}.spline.cost);
    else
        fprintf('[A*] Scene %d: No feasible grid path found（res=%d m）\n',s,cfg.astar_res);
    end
    save(outFile,'data24','data48','astar','done24','done48','doneAstar');
end

%% 
fprintf('\n=====Start rendering=====\n');
fig_uav_results(cfg);
fprintf('\n=====Start making the table=====\n');
table_uav_results(cfg);
fprintf('\n=====All completed=====\n');
end
