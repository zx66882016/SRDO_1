function run_ablations_v6(cfg)


if nargin<1||isempty(cfg), cfg=u_config(); end
outdir=fullfile(cfg.data_dir,'perrun_abl');
if ~exist(outdir,'dir'), mkdir(outdir); end
suite='CEC2017';
funcs=[1 3 5 10 12 15 20 25 28];
stdB=cfg.suites(strcmp({cfg.suites.name},suite)).stdBudget;
budgets=[cfg.fastBudget, stdB];

jobs={};
% ---- 1、Module Ablation ----
mods={
  'Mod_M0','SRDO_NOLR',     struct();
  'Mod_M1','AB_NM',         struct();
  'Mod_M2','AB_CMA',        struct();
  'Mod_M3','AB_CMA_STALL',  struct();
  'Mod_M4','AB_MULTISTART', struct();
  'Mod_M5','AB_RANDTRIG',   struct();
  'Mod_M6','SRDO',          struct()};
for mi=1:size(mods,1)
  o=mods{mi,3}; o.tag=mods{mi,1};
  for bi=1:numel(budgets)
    for f=funcs
      for r=1:cfg.runs
        stem=sprintf('mod_%s_B%d_F%02d_r%02d',mods{mi,1}(5:end),budgets(bi),f,r);
        jobs(end+1,:)={mods{mi,2},suite,f,r,budgets(bi),o,stem}; 
      end
    end
  end
end
% ---- 2、Budget ratio ----
fracs=[0.10 0.18 0.25];
for qi=1:3
  o=struct('localFrac',fracs(qi),'tag',sprintf('Frac_Q%02d',round(fracs(qi)*100)));
  for f=funcs
    for r=1:cfg.runs
      stem=sprintf('frac_Q%02d_B%d_F%02d_r%02d',round(fracs(qi)*100),stdB,f,r);
      jobs(end+1,:)={'AB_AUTO',suite,f,r,stdB,o,stem};
    end
  end
end
% ---- 3、Trigger progress----
trigs=[0.70 0.80 0.90]; ttag={'T70','T80','T90'};
for ti=1:3
  o=struct('localTrigger',trigs(ti),'tag',sprintf('Trig_%s',ttag{ti}));
  for f=funcs
    for r=1:cfg.runs
      stem=sprintf('trig_%s_B%d_F%02d_r%02d',ttag{ti},stdB,f,r);
      jobs(end+1,:)={'AB_AUTO',suite,f,r,stdB,o,stem}; %#ok<AGROW>
    end
  end
end
% ---- 4、CMA initial step size ----
sigs=[0.01 0.02 0.05]; stag={'S01','S02','S05'};
for si=1:3
  o=struct('localSigma',sigs(si),'tag',sprintf('Sig_%s',stag{si}));
  for f=funcs
    for r=1:cfg.runs
      stem=sprintf('sig_%s_B%d_F%02d_r%02d',stag{si},stdB,f,r);
      jobs(end+1,:)={'AB_AUTO',suite,f,r,stdB,o,stem}; %#ok<AGROW>
    end
  end
end

keep=true(size(jobs,1),1);
for j=1:size(jobs,1)
  if exist(fullfile(outdir,[jobs{j,7} '.mat']),'file'), keep(j)=false; end
end
jobs=jobs(keep,:);
NJ=size(jobs,1);
fprintf('[V6 ablation] %d independent runs to be completed (parallelism=%d)\n',NJ,cfg.use_parallel);
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
fprintf('[V6 ablation] completed, failed %d, took %.1f minutes',fail,toc(tStart)/60);
end

function bad=run_one(cfg,outdir,job)
bad=0;
[algo,sname,f,r,budget,o,stem]=deal(job{:});
try
  R=u_dispatch(cfg,sname,f,algo,budget,r,o);
  save(fullfile(outdir,[stem '.mat']),'R','-v7');
catch ME
  bad=1;
  fprintf(2,'[Failed] %s: %s\n',stem,ME.message);
end
end
