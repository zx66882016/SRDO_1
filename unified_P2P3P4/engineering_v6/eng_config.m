function cfg = eng_config(mode)

if nargin < 1, mode = 'full'; end

here = fileparts(mfilename('fullpath'));
cfg.here        = here;
cfg.fw_dir      = fileparts(here);          
cfg.legacy_dir  = fileparts(cfg.fw_dir);    

need = {'SRDO_uni.m','LSHADE_uni.m','IPOP_CMAES_uni.m', ...
        'BudgetObjective.m','algo_style.m','pub_style.m','pub_export.m'};
miss = need(~cellfun(@(f)exist(fullfile(cfg.fw_dir,f),'file')==2, need));
if ~isempty(miss)
    error(['%s was not found in the framework directory.\nPlease move the engineering_v6 folder to the unified_P2P3P4 directory', ...
           '（at the same level as SRDO_uni.m）。'], strjoin(miss,'、'));
end
if exist(fullfile(cfg.legacy_dir,'call_algorithm.m'),'file')~=2
    error('call_algorithm.m not found in the upper directory。');
end
addpath(cfg.fw_dir, cfg.legacy_dir, here);


cfg.seed_base  = 20260000;     
cfg.eng_offset = 900000;      
cfg.N          = 30;          
cfg.maxFEs     = 15000;       
cfg.logPoints  = 200;         
cfg.use_parallel = true;
cfg.n_workers  = max(1, feature('numcores'));


cfg.main_algos   = {'SRDO','SRDO_NOLR','SRDO_HP','LSHADE','IPOP_CMAES'};
cfg.legacy_algos = {'PSO','GA','SAR','DBO','MSO','SFOA','GLO','HBO','BNCO'};
cfg.algos        = [cfg.main_algos, cfg.legacy_algos];
cfg.key_algos    = {'SRDO','SRDO_NOLR','SFOA','GLO','LSHADE','IPOP_CMAES'};


cfg.problems = eng_problems();


cfg.penaltyFactor = 1e6;      
cfg.penaltyMode   = 'quad';    
cfg.feasTol       = 1e-4;      


cfg.data_dir  = fullfile(here,'data','perrun_eng');
cfg.table_dir = fullfile(here,'tables');
cfg.fig_dir   = fullfile(here,'figures');
if ~exist(cfg.data_dir,'dir'),  mkdir(cfg.data_dir);  end
if ~exist(cfg.table_dir,'dir'), mkdir(cfg.table_dir); end
if ~exist(cfg.fig_dir,'dir'),   mkdir(cfg.fig_dir);   end


if strcmp(mode,'smoke')
    cfg.runs       = 2;
    cfg.maxFEs     = 900;
    cfg.use_parallel = false;
    cfg.algos      = {'SRDO','SRDO_NOLR','LSHADE','IPOP_CMAES','PSO','SFOA'};
    cfg.problemIds = [1 3 5];   
    cfg.smoke      = true;
else
    cfg.runs       = 30;
    cfg.problemIds = 1:numel(cfg.problems);
    cfg.smoke      = false;
end
end
