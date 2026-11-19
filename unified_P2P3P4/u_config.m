function cfg = u_config(smoke)


if nargin < 1, smoke = false; end

% -------
here = fileparts(mfilename('fullpath'));
cfg.here        = here;                       
cfg.cec_root    = fileparts(here);           
cfg.data_dir    = fullfile(here, 'data');
cfg.fig_dir     = fullfile(here, 'figures');
cfg.table_dir   = fullfile(here, 'tables');
for d = {cfg.data_dir, cfg.fig_dir, cfg.table_dir}
    if ~exist(d{1}, 'dir'), mkdir(d{1}); end
end

% --------
cfg.seed_base   = 20260000;                   
cfg.runs        = 30;                        
cfg.N           = 30;                         

% ---- CEC----
% CEC2017
cfg.suites(1).name      = 'CEC2017';
cfg.suites(1).funcs     = [1, 3:30];
cfg.suites(1).dim       = 30;
cfg.suites(1).lb        = -100;
cfg.suites(1).ub        = 100;
cfg.suites(1).stdBudget = 300000;
cfg.suites(1).checkpoints = [1000 3000 10000 15000 30000 100000 300000];
cfg.suites(1).fstar     = @(f) 100*f;
% CEC2022：
cfg.suites(2).name      = 'CEC2022';
cfg.suites(2).funcs     = 1:12;
cfg.suites(2).dim       = 20;
cfg.suites(2).lb        = -100;
cfg.suites(2).ub        = 100;
cfg.suites(2).stdBudget = 200000;
cfg.suites(2).checkpoints = [1000 3000 10000 15000 30000 100000 200000];
cec2022_fstar = [300 400 600 800 900 1800 2000 2200 2300 2400 2600 2700];
cfg.suites(2).fstar = @(f) cec2022_fstar(f);

% 
cfg.run_suites   = {'CEC2022'};              
cfg.main_algos   = {'SRDO','SRDO_HP','LSHADE','IPOP_CMAES'};
% 
cfg.include_legacy = true;
cfg.legacy_algos  = {'PSO','GA','SAR','DBO','MSO','SFOA','GLO','HBO','BNCO'};

cfg.fastBudget  = 15000;
cfg.T200Budget  = 6000;                      
cfg.repFunc     = [1 5 10 15 25 28];          
cfg.repFuncNames= {'F1','F5','F10','F15','F25','F28'};
cfg.replaceFunc = [1 5 25];                   

% ----success and ERT ----
cfg.thresholds  = [1e-2 1e-4 1e-6 1e-8];
cfg.successTol  = 1e-8;

% --------
cfg.logPoints   = 200;

% --------
cfg.do_verify   = true;    
cfg.do_main     = true;    
cfg.do_ablation = true;    
cfg.do_check    = true;  
cfg.do_analyze  = true;    
cfg.do_figures  = true;    


cfg.use_parallel = true;
cfg.n_workers    = feature('numcores');   


if smoke
    cfg.runs        = 2;
    cfg.suites(1).funcs     = [1 5];
    cfg.suites(1).checkpoints = [1000 3000 6000];
    cfg.suites(1).stdBudget = 6000;
    cfg.suites(2).funcs     = [1 2];
    cfg.suites(2).checkpoints = [1000 3000 6000];
    cfg.suites(2).stdBudget = 6000;
    cfg.run_suites = {'CEC2017','CEC2022'};
    cfg.main_algos  = {'SRDO','LSHADE','IPOP_CMAES'};
    cfg.include_legacy = true;
    cfg.repFunc     = [1 5];
    cfg.repFuncNames= {'F1','F5'};
    cfg.replaceFunc = [1 5];
    cfg.fastBudget  = 3000;
    cfg.T200Budget  = 1200;
    cfg.logPoints   = 60;
    cfg.smoke = true;
else
    cfg.smoke = false;
end
end
