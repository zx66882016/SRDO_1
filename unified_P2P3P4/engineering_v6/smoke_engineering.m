function smoke_engineering

fprintf('========== 1) Formula fingerprint verification==========\n');
check_engineering;

fprintf('========== 2) smoke run ==========\n');
cfg = eng_config('smoke');
cfg.data_dir  = fullfile(cfg.here,'data','smoke');
cfg.table_dir = fullfile(cfg.here,'data','smoke');
cfg.fig_dir   = fullfile(cfg.here,'data','smoke');
if ~exist(cfg.data_dir,'dir'), mkdir(cfg.data_dir); end
run_engineering_v6(cfg);

fprintf('========== 3) Statistics and Charts ==========\n');
A = analyze_engineering(cfg); 
fig_engineering(cfg);
fprintf('==========The smoke test has been passed, and the formal experiment can be conducted ==========\n');
end
