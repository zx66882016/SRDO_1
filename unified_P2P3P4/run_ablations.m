function run_ablations(cfg)


outdir = fullfile(cfg.data_dir,'perrun_abl');
if ~exist(outdir,'dir'), mkdir(outdir); end
suite='CEC2017';
jobs={};
% ---Parameter Ablation ---
for sc=1:4
    for f=cfg.repFunc
        for r=1:cfg.runs
            o=struct('scheme',sc,'tag',sprintf('Param_S%d_T500',sc));
            jobs(end+1,:)={'param',suite,f,r,cfg.fastBudget,o,sprintf('param_S%d_T500_F%02d_r%02d',sc,f,r)}; %#ok<AGROW>
            o=struct('scheme',sc,'tag',sprintf('Param_S%d_T200',sc));
            jobs(end+1,:)={'param',suite,f,r,cfg.T200Budget,o,sprintf('param_S%d_T200_F%02d_r%02d',sc,f,r)}; %#ok<AGROW>
        end
    end
end
% ---  V0..V4 ---
for v=0:4
    o=struct();
    o.useDiff      = v>=1;
    o.useLens      = v>=2;
    o.useHunger    = v>=3;
    o.useDiversity = v>=4;
    o.usePattern   = v>=4;
    o.tag=sprintf('Oper_V%d',v);
    for f=cfg.repFunc
        for r=1:cfg.runs
            jobs(end+1,:)={'oper',suite,f,r,cfg.fastBudget,o,sprintf('oper_V%d_F%02d_r%02d',v,f,r)}; %#ok<AGROW>
        end
    end
end
% --- R0..R3 ---
for rp=0:3
    o=struct('repl',rp,'tag',sprintf('Repl_R%d',rp));
    for f=cfg.replaceFunc
        for r=1:cfg.runs
            jobs(end+1,:)={'repl',suite,f,r,cfg.fastBudget,o,sprintf('repl_R%d_F%02d_r%02d',rp,f,r)}; 
        end
    end
end

%
keep=true(size(jobs,1),1);
for j=1:size(jobs,1)
    if exist(fullfile(outdir,[jobs{j,7} '.mat']),'file'), keep(j)=false; end
end
jobs=jobs(keep,:);
NJ=size(jobs,1);
fprintf('[Ablation] %d independent runs to be completed (parallelism=%d)\n',NJ,cfg.use_parallel);
tStart=tic; fail=0;
if cfg.use_parallel && NJ>0
    pp=gcp('nocreate');
    if isempty(pp), parpool('local',min(cfg.n_workers,max(1,NJ))); end 
    parfor j=1:NJ
        fail=fail+run_one(cfg,outdir,jobs(j,:)); 
    end
else
    for j=1:NJ
        fail=fail+run_one(cfg,outdir,jobs(j,:));
    end
end
fprintf('[Ablation] completed, failed %d, took %.1f minutes',fail,toc(tStart)/60);
end

function bad=run_one(cfg,outdir,job)
bad=0;
[group,sname,f,r,budget,o,stem]=deal(job{:});
try
    R=u_dispatch(cfg,sname,f,'SRDO_NOLR',budget,r,o);
    R.group=group;
    save(fullfile(outdir,[stem '.mat']),'R','-v7');
catch ME
    bad=1;
    fprintf(2,'[Failed] %s: %s\n',stem,ME.message);
end
end
