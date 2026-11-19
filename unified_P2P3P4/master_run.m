function master_run(smoke)
if nargin<1, smoke=false; end
here=fileparts(mfilename('fullpath'));
addpath(here);
cfg=u_config(smoke);
if smoke, modeStr='smoke test'; else, modeStr='formal experiment'; end
fprintf('==== unified experiment P2/P3/P4，mode：%s ====\n',modeStr);
tAll=tic;
if cfg.do_verify
    ok2=u_verify_optima(cfg);
    if ~ok2
        error('The theoretical optimal point verification has failed. Please check the CEC mex/input_data first and refrain from continuing the experiment。');
    end
end
if cfg.do_main,     run_main_suite(cfg);  end
if cfg.do_ablation, run_ablations(cfg);   end
if cfg.do_check
    ok2=u_check_consistency(cfg);
    if ~ok2, error('P2 cross-table consistency verification failed. Please check and rerun.'); end
end
if cfg.do_analyze
    analyze_main(cfg);
    analyze_ablations(cfg);
end
if cfg.do_figures
    fig_main(cfg);
    fig_ablations(cfg);
end
fprintf('====Completed in full, total time %.1f minutes ====\n',toc(tAll)/60);
fprintf('Table：%s\nFigure：%s\n Sequential data：%s\n',cfg.table_dir,cfg.fig_dir,cfg.data_dir);
end
