function fig_main_box(cfg)
if nargin<1||isempty(cfg), cfg=u_config(); end
keyAlgos={'SRDO','SFOA','GLO','LSHADE','IPOP_CMAES'};
for s=1:numel(cfg.run_suites)
  sname=cfg.run_suites{s};
  si=find(strcmp({cfg.suites.name},sname),1);
  S=cfg.suites(si);
  repf=cfg.repFunc; repf=repf(ismember(repf,S.funcs));
  if strcmp(sname,'CEC2022'), repf=[1 3 6 8 10 12]; repf=repf(ismember(repf,S.funcs)); end
  if numel(repf)>6, repf=repf(1:6); end
  for mode=1:2
    if mode==1, B=cfg.fastBudget; tag='lim'; else, B=S.stdBudget; tag='std'; end
    d=fullfile(cfg.data_dir,sprintf('perrun_main_B%d',B));
    if isempty(dir(fullfile(d,[sname '_*.mat'])))
      fprintf('[Figure] skip box_%s_%s（no B%d data）\n',sname,tag,B); continue;
    end
    E=nan(numel(repf),numel(keyAlgos),cfg.runs);
    for k=1:numel(repf)
      for ai=1:numel(keyAlgos)
        files=dir(fullfile(d,sprintf('%s_F%02d_%s_r*.mat',sname,repf(k),keyAlgos{ai})));
        for r=1:min(numel(files),cfg.runs)
          Q=load(fullfile(d,files(r).name));
          E(k,ai,r)=max(Q.R.finalError,1e-12);
        end
      end
    end
    if all(isnan(E(:))), continue; end
    labs={'SRDO','SFOA','GLO','L-SHADE','IPOP'};
    cols=cellfun(@(a)algo_style(a).color,keyAlgos,'UniformOutput',false);
    fig=figure('Position',[40 40 760 620]);
    tiledlayout(fig,ceil(numel(repf)/3),3,'Padding','compact','TileSpacing','compact');
    for k=1:numel(repf)
      ax=nexttile;
      M=squeeze(E(k,:,:)).';  
      boxplot(ax,max(M,1e-12),'Labels',labs,'Colors',[0.1 0.1 0.1],'Width',0.6);
      hb=findobj(ax,'Tag','Box');
      for i=1:numel(hb)
        idx=numel(hb)-i+1;
        if idx>=1&&idx<=numel(cols)
          patch(get(hb(i),'XData'),get(hb(i),'YData'),cols{idx},...
                'FaceAlpha',0.25,'EdgeColor',cols{idx},'Tag','Box');
        end
      end
      set(ax,'YScale','log'); grid(ax,'on');
      title(ax,sprintf('%s F%d',sname,repf(k)));
      ylabel(ax,'Error f(x)-f*');
      if numel(repf)-k>=3, xlabel(ax,''); end
      pub_style(fig,ax);
      ax.XTickLabelRotation=20;
    end
    pub_export(fig,sprintf('box_%s_%s',sname,tag),cfg); close(fig);
  end
end
end
