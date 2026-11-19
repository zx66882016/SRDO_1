function analyze_runtime(cfg)
if nargin<1||isempty(cfg), cfg=u_config(); end
xlsx=fullfile(cfg.table_dir,'runtime_results.xlsx');
if exist(xlsx,'file'), delete(xlsx); end
for s=1:numel(cfg.run_suites)
  sname=cfg.run_suites{s};
  si=find(strcmp({cfg.suites.name},sname),1);
  S=cfg.suites(si);
  order=cfg.main_algos; if cfg.include_legacy, order=[order,cfg.legacy_algos]; end
  for B=unique([cfg.fastBudget S.stdBudget])
    d=fullfile(cfg.data_dir,sprintf('perrun_main_B%d',B));
    files=dir(fullfile(d,[sname '_*.mat']));
    if isempty(files), continue; end
    Tm=containers.Map();
    for i=1:numel(files)
      tok=regexp(files(i).name,sprintf('%s_F\\d+_(.+)_r\\d+\\.mat',sname),'tokens','once');
      if isempty(tok), continue; end
      algo=tok{1};
      try
        Q=load(fullfile(d,files(i).name));
        if isfield(Q.R,'elapsed') && ~isempty(Q.R.elapsed)
          if ~isKey(Tm,algo), Tm(algo)=[]; end
          Tm(algo)=[Tm(algo); Q.R.elapsed];
        end
      catch
      end
    end
    algos=order(ismember(order,keys(Tm)));
    T=table; T.Algorithm=algos(:);
    med=nan(numel(algos),1); mn=med; q1=med; q3=med; lo=med; hi=med; nr=med;
    for ai=1:numel(algos)
      v=Tm(algos{ai}); med(ai)=median(v); mn(ai)=mean(v);
      q1(ai)=prctile(v,25); q3(ai)=prctile(v,75);
      lo(ai)=min(v); hi(ai)=max(v); nr(ai)=numel(v);
    end
    T.N_runs=nr; T.Median_s=med; T.Mean_s=mn; T.Q1_s=q1; T.Q3_s=q3;
    T.Min_s=lo; T.Max_s=hi;
    T=sortrows(T,'Median_s');
    if B==cfg.fastBudget, bl='limited'; else, bl='standard'; end
    writetable(T,xlsx,'Sheet',sprintf('%s_%s',sname,bl),'WriteMode','overwritesheet');
    fprintf('[machine time] %s @ %d FEs statistics（%d algorithms）\n',sname,B,numel(algos));
  end
end
fprintf('[machine time] already written %s\n',xlsx);
end
