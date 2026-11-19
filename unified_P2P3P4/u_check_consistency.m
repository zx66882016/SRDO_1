function ok = u_check_consistency(cfg)

ok=true; mx=0; n=0;
maindir=fullfile(cfg.data_dir,sprintf('perrun_main_B%d',cfg.fastBudget));
abldir=fullfile(cfg.data_dir,'perrun_abl');
fprintf('\n========== P2 Cross-table consistency verification==========\n');
pairs={{'oper','Oper_V4','oper_V4'},{'param','Param_S1_T500','param_S1_T500'},{'repl','Repl_R0','repl_R0'}};
for fi=1:numel(cfg.repFunc)
    f=cfg.repFunc(fi);
    for r=1:cfg.runs
        mf=fullfile(maindir,sprintf('CEC2017_F%02d_SRDO_NOLR_r%02d.mat',f,r));
        if ~exist(mf,'file'), continue; end
        M=load(mf); a=M.R.finalError;
        for k=1:size(pairs,1)
            fn=[];
            switch pairs{k}{1}
                case 'oper', fn=fullfile(abldir,sprintf('oper_V4_F%02d_r%02d.mat',f,r));
                case 'param',fn=fullfile(abldir,sprintf('param_S1_T500_F%02d_r%02d.mat',f,r));
                case 'repl'
                    if any(f==cfg.replaceFunc)
                        fn=fullfile(abldir,sprintf('repl_R0_F%02d_r%02d.mat',f,r));
                    end
            end
            if isempty(fn)||~exist(fn,'file'), continue; end
            Q=load(fn); d=abs(a-Q.R.finalError);
            mx=max(mx,d); n=n+1;
            if d>0, ok=false; fprintf('  inconsistent F%d r%d vs %s: %.10g vs %.10g\n',f,r,pairs{k}{2},a,Q.R.finalError); end
        end
    end
end
fprintf('Compare %d groups，maximum difference = %.3e\n',n,mx);
if ok, fprintf('=> P2 Pass：The main experiment is fully consistent with the complete SRDO of the three ablation types。\n');
else,  fprintf('=> P2 Not passed：Please check the code version/seed/budget and rerun!\n'); end
end
