function run_all_figs()


here = fileparts(mfilename('fullpath'));
cd(here);
cfg = u_config();

fprintf('\n=== [1/4] Main experimental figure ===\n');
fig_main(cfg);

fprintf('\n=== [2/4] Line diagram of main experimental box ===\n');
fig_main_box(cfg);

fprintf('\n=== [3/4] Old ablation map (operators O0-O4 / parameters S1-S4 / replacements R0-R3) ===\n');
fig_ablations(cfg);

fprintf('\n=== [4/4] V6 Ablation diagram (modules M0-M6 / three sets of sensitivities) ===\n');
fig_v6(cfg);

fprintf('\n=== All drawings have been completed ===\n');
fprintf('Output directory：%s\n', cfg.fig_dir);
fprintf('If a certain image indicates "skip (no data)", it means that the corresponding perrun mat is missing and the experiment for that group needs to be rerun.\n');
end
