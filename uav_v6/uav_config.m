function cfg = uav_config()


here = fileparts(mfilename('fullpath'));
cfg.here       = here;
cfg.data_dir   = fullfile(here,'data');
cfg.fig_dir    = fullfile(here,'figures');
cfg.table_dir  = fullfile(here,'tables');
for d = {cfg.data_dir,cfg.fig_dir,cfg.table_dir}
    if ~exist(d{1},'dir'), mkdir(d{1}); end
end

cfg.path_algos   = fileparts(here);                            
cfg.path_unified = fullfile(fileparts(here),'unified_P2P3P4'); 
addpath(cfg.path_unified);
addpath(cfg.path_algos);


cfg.N      = 30;       
cfg.runs   = 30;       
cfg.fes24  = 12000;   
cfg.fes48  = 24000;   


cfg.scenarios = {'S1_sparse','S2_medium','S3_dense','S4_corridor'};
cfg.hd_scenario = 'S3_dense';   

cfg.algos = {
    'SRDO',       'SRDO';
    'SRDO_NOLR',  'SRDO (w/o LR)';
    'LSHADE',     'L-SHADE';
    'IPOP_CMAES', 'IPOP-CMA-ES';
    'SFOA',       'SFOA';
    'GLO',        'GLO';
    'PSO',        'PSO';
    'GA',         'GA';
    'SAR',        'SAR';
    'DBO',        'DBO'};



cfg.astar_res = 5;


cfg.space      = 100;    
cfg.n_sample   = 50;     
cfg.safe_dist  = 5;      
cfg.w_len      = 1.0;    
cfg.w_obs      = 0.01;   
cfg.w_smooth   = 0.5;  
cfg.w_mono     = 0.1;   
cfg.pen_hard   = 1e6;    
cfg.pen_soft   = 1e3;    
end
