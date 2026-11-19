function run_main_suite(cfg)


algos = cfg.main_algos;
if cfg.include_legacy, algos = [algos, cfg.legacy_algos]; end


jobs={};
for s=1:numel(cfg.run_suites)
    sname=cfg.run_suites{s};
    si=find(strcmp({cfg.suites.name},sname),1);
    S=cfg.suites(si);
    budgets=unique([cfg.fastBudget, S.stdBudget]);
    for bi=1:numel(budgets)
        d=fullfile(cfg.data_dir,sprintf('perrun_main_B%d',budgets(bi)));
        if ~exist(d,'dir'), mkdir(d); end
        for a=1:numel(algos)
            for f=S.funcs
                for r=1:cfg.runs
                    fn=fullfile(d,sprintf('%s_F%02d_%s_r%02d.mat',sname,f,algos{a},r));
                    if ~exist(fn,'file')
                        jobs(end+1,:)={sname,budgets(bi),f,algos{a},r,fn}; 
                    end
                end
            end
        end
    end
end
NJ=size(jobs,1);
fprintf('[Main experiment] %d independent runs to be completed (parallelism=%d)\n',NJ,cfg.use_parallel);
tStart=tic;
fail=0;
if cfg.use_parallel && NJ>0
    pp=gcp('nocreate');
    if isempty(pp), parpool('local',min(cfg.n_workers,max(1,NJ))); end 
    parfor j=1:NJ
        fail=fail+run_one(cfg,jobs(j,:)); 
    end
else
    for j=1:NJ
        fail=fail+run_one(cfg,jobs(j,:));
    end
end
fprintf('[Main experiment] completed, failed %d, took %.1f minutes\n',fail,toc(tStart)/60);
end

function bad=run_one(cfg,job)
bad=0;
[sname,budget,f,algo,r,fn]=deal(job{:});
try
    R=u_dispatch(cfg,sname,f,algo,budget,r);
    save(fn,'R','-v7');
    fprintf('  %s B%d F%d %s r%d FE=%d err=%.3g (%.1fs)\n',...
        sname,budget,f,algo,r,R.nFEs,R.finalError,R.elapsed);
catch ME
    bad=1;
    fprintf(2,'[Failed] %s B%d F%d %s r%d: %s\n',sname,budget,f,algo,r,ME.message);
end
end
