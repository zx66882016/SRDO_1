function run_engineering_v6(cfg)

if nargin < 1, cfg = eng_config(); end
if ischar(cfg) || isstring(cfg), cfg = eng_config(cfg); end

algos = cfg.algos;
pids  = cfg.problemIds;
runs  = cfg.runs;

jobs = {};
for ip = 1:numel(pids)
    pid = pids(ip);
    for a = 1:numel(algos)
        for r = 1:runs
            fn = fullfile(cfg.data_dir, sprintf('%s_%s_r%02d.mat', ...
                cfg.problems(pid).tag, algos{a}, r));
            if ~exist(fn,'file')
                jobs(end+1,:) = {pid, algos{a}, r, fn}; 
            end
        end
    end
end
NJ = size(jobs,1);
fprintf('[Engineering Experiment] %d independent runs to be completed (parallelism=%d, budget %d FEs/run)\n', ...
    NJ, cfg.use_parallel, cfg.maxFEs);
tStart = tic;
fail = 0;
if cfg.use_parallel && NJ > 0
    pp = gcp('nocreate');
    if isempty(pp), parpool('local', min(cfg.n_workers, max(1,NJ))); end 
    parfor j = 1:NJ
        fail = fail + run_one(cfg, jobs(j,:)); 
    end
else
    for j = 1:NJ
        fail = fail + run_one(cfg, jobs(j,:));
    end
end
fprintf('[Engineering experiment] completed, failed %d, took %.1f minutes', fail, toc(tStart)/60);
end

function bad = run_one(cfg, job)
bad = 0;
[pid, algo, r, fn] = deal(job{:});
try
    R = eng_dispatch(cfg, pid, algo, r);
    save(fn,'R','-v7');
    fprintf('  P%d %s %s r%02d FE=%d f=%.6g feas=%d (%.1fs)\n', ...
        pid, R.tag, algo, r, R.nFEs, R.finalRawAll, R.feasible, R.elapsed);
catch ME
    bad = 1;
    fprintf(2,'[failed] P%d %s r%02d: %s\n', pid, algo, r, ME.message);
end
end
