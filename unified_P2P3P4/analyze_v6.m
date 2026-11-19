function AV = analyze_v6(cfg)
if nargin<1||isempty(cfg), cfg=u_config(); end
d=fullfile(cfg.data_dir,'perrun_abl');
xlsx=fullfile(cfg.table_dir,'v6_ablation_results.xlsx');
if exist(xlsx,'file'), delete(xlsx); end
funcs=[1 3 5 10 12 15 20 25 28];
stdB=cfg.suites(strcmp({cfg.suites.name},'CEC2017')).stdBudget;
limB=cfg.fastBudget;
AV=struct();

% ----Module Ablation----
mTags={'M0','M1','M2','M3','M4','M5','M6'};
for mode=1:2
  if mode==1, B=limB; sheet='mod_limited'; else, B=stdB; sheet='mod_standard'; end
  E=nan(numel(funcs),7,cfg.runs);
  for fi=1:numel(funcs)
    for mi=1:7
      for r=1:cfg.runs
        fn=fullfile(d,sprintf('mod_%s_B%d_F%02d_r%02d.mat',mTags{mi},B,funcs(fi),r));
        if exist(fn,'file'), S=load(fn); E(fi,mi,r)=S.R.finalError; end
      end
    end
  end
  Mean=squeeze(mean(E,3,'omitnan')); Std=squeeze(std(E,0,3,'omitnan'));
  Min=squeeze(min(E,[],3));
  rk=zeros(size(Mean));
  for fi=1:numel(funcs), rk(fi,:)=tr(Mean(fi,:)); end
  T=table; T.Func=arrayfun(@(f)sprintf('F%d',f),funcs,'UniformOutput',false).';
  for mi=1:7
    T.([mTags{mi} '_mean'])=Mean(:,mi); T.([mTags{mi} '_std'])=Std(:,mi);
    T.([mTags{mi} '_min'])=Min(:,mi);   T.([mTags{mi} '_rank'])=rk(:,mi);
  end
  % Wilcoxon
  for base=[1 6]
    others=setdiff(1:7,base); P=nan(numel(funcs),7); MK=strings(numel(funcs),7);
    for fi=1:numel(funcs)
      ps=arrayfun(@(mi) sr(squeeze(E(fi,base,:)),squeeze(E(fi,mi,:))),others);
      [~,~,pa]=holm(ps);
      for k=1:numel(others)
        mi=others(k); P(fi,mi)=pa(k);
        MK(fi,mi)=mkmark(pa(k),Mean(fi,base),Mean(fi,mi));
      end
    end
    for mi=1:7
      T.(sprintf('p_%s_vs_%s',mTags{mi},mTags{base}))=P(:,mi);
      T.(sprintf('mark_%s_vs_%s',mTags{mi},mTags{base}))=cellstr(MK(:,mi));
    end
  end
  writetable(T,xlsx,'Sheet',sheet,'WriteMode','overwritesheet');
  AV.(sheet)=struct('funcs',funcs,'tags',{mTags},'E',E,'mean',Mean,'std',Std,...
                    'min',Min,'rank',rk,'avgRank',mean(rk,1),'budget',B);
end

% ----Three types of parameter sensitivities (standard budget)----
grp=repmat(struct('sheet','','prefix','','labs',{{}},'vals',[],'mid',1),1,3);
grp(1)=struct('sheet','frac_standard','prefix','frac_Q','labs',{{'Q10','Q18','Q25'}},...
              'vals',[10 18 25],'mid',2);
grp(2)=struct('sheet','trig_standard','prefix','trig_T','labs',{{'T70','T80','T90'}},...
              'vals',[70 80 90],'mid',2);
grp(3)=struct('sheet','sig_standard','prefix','sig_S','labs',{{'S01','S02','S05'}},...
              'vals',[1 2 5],'mid',2);
for g=1:3
  n=numel(grp(g).labs); E=nan(numel(funcs),n,cfg.runs);
  for fi=1:numel(funcs)
    for k=1:n
      for r=1:cfg.runs
        fn=fullfile(d,sprintf('%s%02d_B%d_F%02d_r%02d.mat',grp(g).prefix,grp(g).vals(k),stdB,funcs(fi),r));
        if exist(fn,'file'), S=load(fn); E(fi,k,r)=S.R.finalError; end
      end
    end
  end
  Mean=squeeze(mean(E,3,'omitnan')); Std=squeeze(std(E,0,3,'omitnan'));
  Min=squeeze(min(E,[],3));
  rk=zeros(size(Mean));
  for fi=1:numel(funcs), rk(fi,:)=tr(Mean(fi,:)); end
  T=table; T.Func=arrayfun(@(f)sprintf('F%d',f),funcs,'UniformOutput',false).';
  for k=1:n
    T.([grp(g).labs{k} '_mean'])=Mean(:,k); T.([grp(g).labs{k} '_std'])=Std(:,k);
    T.([grp(g).labs{k} '_min'])=Min(:,k);    T.([grp(g).labs{k} '_rank'])=rk(:,k);
  end
  others=setdiff(1:n,grp(g).mid); P=nan(numel(funcs),n); MK=strings(numel(funcs),n);
  for fi=1:numel(funcs)
    ps=arrayfun(@(k) sr(squeeze(E(fi,grp(g).mid,:)),squeeze(E(fi,k,:))),others);
    [~,~,pa]=holm(ps);
    for ii=1:numel(others)
      k=others(ii); P(fi,k)=pa(ii);
      MK(fi,k)=mkmark(pa(ii),Mean(fi,grp(g).mid),Mean(fi,k));
    end
  end
  for k=1:n
    T.(sprintf('p_%s_vs_%s',grp(g).labs{k},grp(g).labs{grp(g).mid}))=P(:,k);
    T.(sprintf('mark_%s_vs_%s',grp(g).labs{k},grp(g).labs{grp(g).mid}))=cellstr(MK(:,k));
  end
  writetable(T,xlsx,'Sheet',grp(g).sheet,'WriteMode','overwritesheet');
  AV.(grp(g).sheet)=struct('funcs',funcs,'labs',{grp(g).labs},'E',E,'mean',Mean,...
                           'std',Std,'min',Min,'rank',rk,'avgRank',mean(rk,1));
end

save(fullfile(cfg.data_dir,'analysis_v6.mat'),'AV','-v7');
fprintf('[V6 statistics] already written %s\n',xlsx);
fprintf('Marking description: + indicates significantly better control; - indicates significantly better variant; = indicates no significant difference (after Holm correction, p<0.05)\n');
end

function r=tr(x)
if exist('tiedrank','file')==2, r=tiedrank(x); else
  [~,ix]=sort(x); r=zeros(size(x)); r(ix)=1:numel(x);
  u=unique(x); for i=1:numel(u), r(x==u(i))=mean(r(x==u(i))); end
end
end
function p=sr(a,b)
if exist('signrank','file')==2 && numel(a)==numel(b) && any(~isnan(a)) && any(~isnan(b))
  p=signrank(a,b);
else
  p=NaN;
end
end
function [h,padj,praw]=holm(p)
m=numel(p); praw=p(:); [ps,ix]=sort(praw); adj=ps.*(m-(1:m).'); adj=min(adj,1);
for i=2:m, adj(i)=max(adj(i),adj(i-1)); end
padj=nan(m,1); padj(ix)=adj; h=padj<0.05;
end
function mk=mkmark(p,mc,mv)
if isnan(p), mk='n/a'; elseif p>=0.05, mk='=';
elseif mc<mv, mk='+'; else, mk='-'; end
end
