function reset_results()


cfg=u_config(false);
dirs={fullfile(cfg.data_dir,'perrun_main'),...
      fullfile(cfg.data_dir,'perrun_abl')};
d=dir(fullfile(cfg.data_dir,'perrun_main_B*'));
for i=1:numel(d), dirs{end+1}=fullfile(cfg.data_dir,d(i).name); end 
for i=1:numel(dirs)
    if exist(dirs{i},'dir'), rmdir(dirs{i},'s'); fprintf('Deleted %s\n',dirs{i}); end
end
for f={fullfile(cfg.data_dir,'analysis_main.mat'),fullfile(cfg.data_dir,'analysis_abl.mat')}
    if exist(f{1},'file'), delete(f{1}); fprintf('Deleted %s\n',f{1}); end
end
tt=dir(fullfile(cfg.table_dir,'*'));
for i=1:numel(tt), if ~tt(i).isdir, delete(fullfile(cfg.table_dir,tt(i).name)); end; end
ff=dir(fullfile(cfg.fig_dir,'*'));
for i=1:numel(ff), if ~ff(i).isdir, delete(fullfile(cfg.fig_dir,ff(i).name)); end; end
fprintf('Cleanup completed, you can now officially start running master_run。\n');
end
