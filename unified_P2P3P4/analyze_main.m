function A = analyze_main(cfg)
xlsx=fullfile(cfg.table_dir,'main_results.xlsx');
if exist(xlsx,'file'), delete(xlsx); end
A=struct();
for s=1:numel(cfg.run_suites)
    sname=cfg.run_suites{s};
    si=find(strcmp({cfg.suites.name},sname),1);
    S=cfg.suites(si);
    funcs=S.funcs; stdB=S.stdBudget;
    budgets=[cfg.fastBudget, stdB]; blabels={'limited','standard'};
    % 
    tmp=dir(fullfile(cfg.data_dir,sprintf('perrun_main_B%d',stdB),'*.mat'));
    if isempty(tmp), warning('Not found %s Main experimental results，skip',sname); continue; end
    algos={};
    for i=1:numel(tmp)
        tok=regexp(tmp(i).name,sprintf('%s_F\\d+_(.+)_r\\d+\\.mat',sname),'tokens','once');
        if ~isempty(tok), algos{end+1}=tok{1}; end 
    end
    order=cfg.main_algos; if cfg.include_legacy, order=[order,cfg.legacy_algos]; end
    algos=unique(algos,'stable');
    algos=order(ismember(order,algos));
    for bi=1:2
        B=budgets(bi);
        R=loadR(fullfile(cfg.data_dir,sprintf('perrun_main_B%d',B)));
        E=nan(numel(funcs),numel(algos),cfg.runs);
        for fi=1:numel(funcs)
            for ai=1:numel(algos)
                key=sprintf('%s|%d|%s',sname,funcs(fi),algos{ai});
                if ~isfield(R,mangle(key)), continue; end
                rr=R.(mangle(key)).runs;
                for k=1:min(numel(rr),cfg.runs)
                    idx=find(rr{k}.cpFEs>=B,1);
                    E(fi,ai,k)=rr{k}.cpError(idx);
                end
            end
        end
        Mean=squeeze(mean(E,3,'omitnan')); Std=squeeze(std(E,0,3,'omitnan'));
        SEM=Std/sqrt(cfg.runs);
        meanE=squeeze(mean(E,3,'omitnan'));
        rk=zeros(size(meanE));
        for fi=1:numel(funcs), rk(fi,:)=tr(meanE(fi,:)); end
        avgRank=mean(rk,1);

        T=table; T.Func=arrayfun(@(f)sprintf('F%d',f),funcs,'UniformOutput',false).';
        for ai=1:numel(algos)
            T.([algos{ai} '_mean'])=Mean(:,ai);
            T.([algos{ai} '_std'])=Std(:,ai);
            T.([algos{ai} '_SEM'])=SEM(:,ai);
            T.([algos{ai} '_rank'])=rk(:,ai);
        end
        sheet=[sname '_' blabels{bi}];
        writetable(T,xlsx,'Sheet',sheet,'WriteMode','overwritesheet');
        RankT=table(algos(:),avgRank(:),'VariableNames',{'Algorithm','AvgFriedmanRank'});
        writetable(RankT,xlsx,'Sheet',[sheet '_rank'],'WriteMode','overwritesheet');

        base=find(strcmp(algos,'SRDO'),1);
        pmat=nan(numel(funcs),numel(algos)); hmat=strings(numel(funcs),numel(algos));
        a12=nan(numel(funcs),numel(algos));
        for fi=1:numel(funcs)
            ps=[]; pairs=[];
            for ai=1:numel(algos)
                if ai==base, continue; end
                x=squeeze(E(fi,base,:)); y=squeeze(E(fi,ai,:));
                a12(fi,ai)=A12(x,y);
                if exist('signrank','file')==2 && numel(x)==numel(y), p=signrank(x,y); else, p=NaN; end
                pmat(fi,ai)=p; ps(end+1)=p; pairs(end+1)=ai; 
            end
            [~,~,padj]=holm(ps);
            for kk=1:numel(pairs)
                ai=pairs(kk);
                hmat(fi,ai)=sigmark(padj(kk),squeeze(E(fi,base,:)),squeeze(E(fi,ai,:)));
                pmat(fi,ai)=padj(kk);
            end
            hmat(fi,base)='-';
        end
        WT=table; WT.Func=T.Func;
        for ai=1:numel(algos)
            WT.([algos{ai} '_p'])=pmat(:,ai);
            WT.([algos{ai} '_mark'])=cellstr(hmat(:,ai));
            WT.([algos{ai} '_A12'])=a12(:,ai);
        end
        writetable(WT,xlsx,'Sheet',[sheet '_wilcoxon'],'WriteMode','overwritesheet');
        A.([sname '_' blabels{bi}])=struct('funcs',funcs,'algos',{algos},'E',E,...
            'mean',Mean,'std',Std,'sem',SEM,'rank',rk,'avgRank',avgRank,...
            'p',pmat,'mark',hmat,'a12',a12,'budget',B);
    end

    % ---- ERT / success rate----
    R=loadR(fullfile(cfg.data_dir,sprintf('perrun_main_B%d',stdB)));
    th=cfg.thresholds; SR=zeros(numel(algos),numel(th)); ERT=nan(numel(algos),numel(th));
    for ai=1:numel(algos)
        for ti=1:numel(th)
            hit=[]; fehit=[];
            for fi=1:numel(funcs)
                key=sprintf('%s|%d|%s',sname,funcs(fi),algos{ai});
                if ~isfield(R,mangle(key)), continue; end
                rr=R.(mangle(key)).runs;
                for k=1:min(numel(rr),cfg.runs)
                    lx=rr{k}.logFEs; ly=rr{k}.logError;
                    h=find(ly<=th(ti),1);
                    if ~isempty(h), hit(end+1,1)=1; fehit(end+1,1)=lx(h); 
                    else, hit(end+1,1)=0; fehit(end+1,1)=NaN; end         
                end
            end
            ns=sum(hit); SR(ai,ti)=ns/max(1,numel(hit));
            if ns>0, fehit(isnan(fehit))=stdB; ERT(ai,ti)=sum(fehit)/ns; end
        end
    end
    ET=table(algos(:),SR,ERT,'VariableNames',{'Algorithm','SuccessRate','ERT'});
    writetable(ET,xlsx,'Sheet',[sname '_ERT'],'WriteMode','overwritesheet');
    A.([sname '_ERT'])=struct('algos',{algos},'th',th,'SR',SR,'ERT',ERT);
    save(fullfile(cfg.data_dir,'analysis_main.mat'),'A','-v7');
end
fprintf('[Statistics] Already written %s\n',xlsx);
end
% ---------------- tool ----------------
function R=loadR(d)
R=struct();
files=dir(fullfile(d,'*.mat'));
for i=1:numel(files)
    S=load(fullfile(d,files(i).name)); q=S.R;
    key=sprintf('%s|%d|%s',q.suite,q.func,q.algo);
    if ~isfield(R,mangle(key))
        R.(mangle(key))=struct('runs',{{}});
    end
    R.(mangle(key)).runs{end+1}=q;
end
end
function s=mangle(k), s=['k' char(regexprep(k,'[^a-zA-Z0-9]','_'))]; end
function r=tr(x)
if exist('tiedrank','file')==2, r=tiedrank(x); else
    [~,ix]=sort(x); r=zeros(size(x)); r(ix)=1:numel(x);
    u=unique(x); for i=1:numel(u), r(x==u(i))=mean(r(x==u(i))); end
end
end
function a=A12(x,y)
x=x(:); y=y(:); m=numel(x); n=numel(y);
a=(sum(x>y)+0.5*sum(x==y))/(m*n);
end
function [h,padj,praw]=holm(p)
m=numel(p); praw=p(:); [ps,ix]=sort(praw); adj=ps.*(m-(1:m).'); adj=min(adj,1);
for i=2:m, adj(i)=max(adj(i),adj(i-1)); end
padj=nan(m,1); padj(ix)=adj; h=padj<0.05;
end
function mk=sigmark(p,x,y)
if isnan(p), mk='n/a'; elseif p>=0.05, mk='=';
elseif mean(x)<mean(y), mk='+'; else, mk='-'; end
end