function fig_v6(cfg)
if nargin<1||isempty(cfg), cfg=u_config(); end
d=fullfile(cfg.data_dir,'perrun_abl');
funcs=[1 3 5 10 12 15 20 25 28];
stdB=cfg.suites(strcmp({cfg.suites.name},'CEC2017')).stdBudget;
limB=cfg.fastBudget;

mAlgo={'M0','M1','M2','M3','M4','M5','M6'};
for mode=1:2
  if mode==1, B=limB; tag='lim'; else, B=stdB; tag='std'; end
  make_group(cfg,d,funcs,B,sprintf('mod_%s_B%d',mAlgo{1},B),mAlgo,'M',...
             sprintf('abl_v6mod_conv_%s',tag),'conv');
  make_group(cfg,d,funcs,B,sprintf('mod_%s_B%d',mAlgo{1},B),mAlgo,'M',...
             sprintf('abl_v6mod_box_%s',tag),'box');
end

make_group(cfg,d,funcs,stdB,sprintf('frac_Q10_B%d',stdB),{'Q10','Q18','Q25'},'frac','abl_v6_frac_conv','conv');
make_group(cfg,d,funcs,stdB,sprintf('frac_Q10_B%d',stdB),{'Q10','Q18','Q25'},'frac','abl_v6_frac_box','box');
make_group(cfg,d,funcs,stdB,sprintf('trig_T70_B%d',stdB),{'T70','T80','T90'},'trig','abl_v6_trig_conv','conv');
make_group(cfg,d,funcs,stdB,sprintf('trig_T70_B%d',stdB),{'T70','T80','T90'},'trig','abl_v6_trig_box','box');
make_group(cfg,d,funcs,stdB,sprintf('sig_S01_B%d',stdB),{'S01','S02','S05'},'sig','abl_v6_sig_conv','conv');
make_group(cfg,d,funcs,stdB,sprintf('sig_S01_B%d',stdB),{'S01','S02','S05'},'sig','abl_v6_sig_box','box');
end

function make_group(cfg,d,funcs,B,probe,labs,kind,figname,kindPlot)
[xx,~]=loadlogs(d,kind,labs{1},funcs(1),B);
if isempty(xx), fprintf('[Figure] skip %s（no B%d data）\n',figname,B); return; end
n=numel(labs);
fig=figure('Position',[40 40 760 640]);
tiledlayout(fig,3,3,'Padding','compact','TileSpacing','compact');
for k=1:numel(funcs)
  f=funcs(k); ax=nexttile; hold(ax,'on');
  if strcmp(kindPlot,'conv')
    anyLine=false;
    for i=1:n
      [x,Y]=loadlogs(d,kind,labs{i},f,B);
      if isempty(x), continue; end
      st=v6_style(kind,labs{i});
      semilogy(ax,x,max(median(Y,1),1e-12),'Color',st.color,'LineStyle',st.ls,...
               'LineWidth',st.lw,'Marker',st.mk,'MarkerSize',3.5,...
               'MarkerIndices',round(linspace(1,numel(x),6)),'DisplayName',st.disp);
      anyLine=true;
    end
    if anyLine
      grid(ax,'on'); xlabel(ax,'Function evaluations'); ylabel(ax,'Error f(x)-f*');
      xlim(ax,[1 B]);
    end
  else
    if exist('boxplot','file')~=2
      text(ax,0.5,0.5,'Statistics Toolbox missing（no boxplot）','HorizontalAlignment','center');
    else
      M=nan(cfg.runs,n);
      for i=1:n
        Y=loadend(d,kind,labs{i},f,B);
        if ~isempty(Y), M(1:numel(Y),i)=max(Y(:),1e-12); end
      end
      if all(isnan(M(:))), continue; end
      cols=arrayfun(@(i)v6_style(kind,labs{i}).color,1:n,'UniformOutput',false);
      boxplot(ax,max(M,1e-12),'Labels',labs,'Colors',[0.1 0.1 0.1],'Width',0.6);
      
      hb=findobj(ax,'Tag','Box');
      for i=1:numel(hb)
        idx=n-i+1;  
        if idx>=1&&idx<=n, patch(get(hb(i),'XData'),get(hb(i),'YData'),cols{idx},...
            'FaceAlpha',0.25,'EdgeColor',cols{idx},'Tag','Box'); end
      end
      set(ax,'YScale','log'); grid(ax,'on');
    end
  end
  title(ax,sprintf('F%d',f)); pub_style(fig,ax);
  if k==1
    if strcmp(kindPlot,'conv'), legend(ax,'Location','northeast');
    else, legend(ax,'off'); end
  end
end
pub_export(fig,figname,cfg); close(fig);
end

function [x,Y]=loadlogs(d,kind,lab,f,B)
x=[]; Y=[];
switch kind
  case 'M',    pat=sprintf('mod_%s_B%d_F%02d_r*.mat',lab,B,f);
  case 'frac', pat=sprintf('frac_%s_B%d_F%02d_r*.mat',lab,B,f);
  case 'trig', pat=sprintf('trig_%s_B%d_F%02d_r*.mat',lab,B,f);
  case 'sig',  pat=sprintf('sig_%s_B%d_F%02d_r*.mat',lab,B,f);
end
files=dir(fullfile(d,pat));
for i=1:numel(files)
  S=load(fullfile(d,files(i).name)); q=S.R;
  if isempty(x), x=q.logFEs; end
  Y(end+1,:)=q.logError; 
end
end

function Y=loadend(d,kind,lab,f,B)
Y=[];
switch kind
  case 'M',    pat=sprintf('mod_%s_B%d_F%02d_r*.mat',lab,B,f);
  case 'frac', pat=sprintf('frac_%s_B%d_F%02d_r*.mat',lab,B,f);
  case 'trig', pat=sprintf('trig_%s_B%d_F%02d_r*.mat',lab,B,f);
  case 'sig',  pat=sprintf('sig_%s_B%d_F%02d_r*.mat',lab,B,f);
end
files=dir(fullfile(d,pat));
for i=1:numel(files)
  S=load(fullfile(d,files(i).name));
  Y(end+1,1)=S.R.finalError; 
end
end

function st=v6_style(kind,lab)

modC=[0.35 0.35 0.35; 0.00 0.62 0.45; 0.20 0.55 0.85; 0.55 0.25 0.70; ...
      0.84 0.37 0.00; 0.50 0.50 0.50; 0.80 0.15 0.15];
modL={':','--','--',':','-.',':','-'};
modW=[1.2 1.3 1.3 1.3 1.3 1.3 2.2];
modM={'none','none','none','none','none','none','pentagram'};
modD={'SRDO (w/o LR)','NM-only (M1)','CMA-only (M2)','Stall-trig (M3)',...
      'Multi-start (M4)','Random-trig (M5)','SRDO'};
senC=[0.35 0.70 0.90; 0.80 0.15 0.15; 0.20 0.30 0.65];
senL={':','-','--'}; senM={'none','pentagram','none'};
switch kind
  case 'M'
    i=find(strcmp({'M0','M1','M2','M3','M4','M5','M6'},lab),1);
    st=struct('color',modC(i,:),'ls',modL{i},'lw',modW(i),'mk',modM{i},'disp',modD{i});
  otherwise
    switch lab
      case {'Q10','T70','S01'}, i=1;
      case {'Q18','T80','S02'}, i=2;
      otherwise, i=3;
    end
    st=struct('color',senC(i,:),'ls',senL{i},'lw',1.4+0.8*(i==2),...
              'mk',senM{i},'disp',lab);
end
end
