function A = analyze_engineering(cfg)


if nargin < 1, cfg = eng_config(); end
algos = cfg.algos; pids = cfg.problemIds; R = cfg.runs;
nP = numel(pids); nA = numel(algos);

% --------------------------------
F = nan(nP,nA,R);       
Fall = nan(nP,nA,R);     
FEAS = false(nP,nA,R); VIOL = nan(nP,nA,R); TM = nan(nP,nA,R);
bestDesign = cell(nP,nA);
logGrid = []; medCurve = cell(nP,nA); allCurve = cell(nP,nA);

for ip = 1:nP
    pid = pids(ip); prob = cfg.problems(pid);
    for ai = 1:nA
        vals = nan(R,1); fallv = nan(R,1); fs = false(R,1);
        vv = nan(R,1); tt = nan(R,1); curves = [];
        bestf = inf; bestR = [];
        for r = 1:R
            fn = fullfile(cfg.data_dir, sprintf('%s_%s_r%02d.mat',prob.tag,algos{ai},r));
            if ~exist(fn,'file'), continue; end
            S = load(fn); q = S.R;
            fallv(r) = q.finalRawAll; fs(r) = q.feasible;
            vv(r) = q.viol; tt(r) = q.elapsed;
            if q.feasible
                vals(r) = q.finalRaw;
                if q.finalRaw < bestf, bestf = q.finalRaw; bestR = q; end
            end
            if isempty(logGrid), logGrid = q.logFEs; end
            c = q.logBestFeas(:).';
            if isempty(curves), curves = nan(R,numel(c)); end
            curves(r,:) = c(1:size(curves,2));
        end
        F(ip,ai,:) = vals; Fall(ip,ai,:) = fallv; FEAS(ip,ai,:)=fs;
        VIOL(ip,ai,:)=vv; TM(ip,ai,:)=tt;
        if ~isempty(bestR), bestDesign{ip,ai} = bestR; end
        allCurve{ip,ai} = curves;
        medCurve{ip,ai} = median(curves,1,'omitnan');
    end
end

% -------------------------------
Mean = squeeze(mean(F,3,'omitnan'));
Std  = squeeze(std(F,0,3,'omitnan'));
SEM  = Std/sqrt(R);
Med  = squeeze(median(F,3,'omitnan'));
Best = squeeze(min(F,[],3));
Worst= squeeze(max(F,[],3));
FeasRate = squeeze(mean(FEAS,3));
MeanViol= squeeze(mean(VIOL,3,'omitnan'));
MeanTime= squeeze(mean(TM,3,'omitnan'));
rk = zeros(nP,nA);
for ip = 1:nP
    v = Mean(ip,:);
    if all(isnan(v)), rk(ip,:)=nA; else
        vv=v; vv(isnan(vv))=max(v(~isnan(v)))+1e6;
        rk(ip,:)=tr(vv);
    end
end
avgRank = mean(rk,1);

% --------------------------------
rowsT = cell(nP*nA,1); rt=0;
for ip=1:nP
    prob=cfg.problems(pids(ip));
    for ai=1:nA
        row=struct();
        row.Problem={prob.tag}; row.Algorithm={algos{ai}};
        row.Best=Best(ip,ai); row.Worst=Worst(ip,ai); row.Mean=Mean(ip,ai);
        row.Std=Std(ip,ai); row.SEM=SEM(ip,ai); row.Median=Med(ip,ai);
        row.FeasRate=FeasRate(ip,ai); row.MeanViol=MeanViol(ip,ai);
        row.MeanTime_s=MeanTime(ip,ai); row.Rank=rk(ip,ai);
        rt=rt+1; rowsT{rt}=struct2table(row);
    end
end
T=vertcat(rowsT{:});

% ---------------- Sheet: Friedman ----------------
TR = table(algos(:), avgRank(:), 'VariableNames', {'Algorithm','AvgFriedmanRank'});
TR = sortrows(TR,'AvgFriedmanRank');

% ---------------- Wilcoxon+ A12 ----------------
[MK6,P6,A12_6] = pairtests(F,'SRDO',algos);
[MK0,P0,A12_0] = pairtests(F,'SRDO_NOLR',algos);
pnames = arrayfun(@(ip)cfg.problems(pids(ip)).tag, 1:nP, 'UniformOutput', false);
TW6 = wtbl(pnames,algos,MK6,P6,A12_6,'SRDO');
TW0 = wtbl(pnames,algos,MK0,P0,A12_0,'SRDO_NOLR');

% ---------------- Sheet: Feasibility ----------------
TF = array2table(FeasRate,'RowNames',pnames,'VariableNames',algos);
TF = rows2vars(TF); TF.Properties.VariableNames{1}='Problem';

% ---------------- Sheet: BestDesign_SRDO ----------------
maxDim=max(arrayfun(@(ip)cfg.problems(pids(ip)).dim,1:nP));
maxCon=max(arrayfun(@(ip)numel(cfg.problems(pids(ip)).cScale),1:nP));
rowsB=cell(nP,1);
sri = find(strcmp(algos,'SRDO'),1);
for ip=1:nP
    prob=cfg.problems(pids(ip)); q=bestDesign{ip,sri};
    row=struct(); row.Problem={prob.tag};
    xv=nan(1,maxDim); gv=nan(1,maxCon);
    if isempty(q)
        row.f=NaN; row.f_ref=prob.fref; row.gap_pct=NaN; row.feasible=false;
    else
        row.f=q.finalRaw; row.f_ref=prob.fref;
        row.gap_pct=(q.finalRaw-prob.fref)/abs(prob.fref)*100;
        row.feasible=q.feasible;
        xv(1:prob.dim)=q.bestXraw; gv(1:numel(q.g))=q.g;
    end
    for j=1:maxDim, row.(sprintf('x%d',j))=xv(j); end
    for j=1:maxCon, row.(sprintf('g%d',j))=gv(j); end
    rowsB{ip}=struct2table(row);
end
TB=vertcat(rowsB{:});

xlsx = fullfile(cfg.table_dir,'engineering_results.xlsx');
if exist(xlsx,'file'), delete(xlsx); end
writetable(T,  xlsx,'Sheet','Stats','WriteMode','overwritesheet');
writetable(TR, xlsx,'Sheet','Friedman','WriteMode','overwritesheet');
writetable(TW6, xlsx,'Sheet','Wilcoxon_SRDO','WriteMode','overwritesheet');
writetable(TW0, xlsx,'Sheet','Wilcoxon_NOLR','WriteMode','overwritesheet');
writetable(TF,  xlsx,'Sheet','Feasibility','WriteMode','overwritesheet');
writetable(TB,  xlsx,'Sheet','BestDesign_SRDO','WriteMode','overwritesheet');

A = struct('algos',{algos},'pids',pids,'pnames',{pnames}, ...
    'F',F,'Fall',Fall,'FEAS',FEAS,'Mean',Mean,'Std',Std,'SEM',SEM, ...
    'Med',Med,'Best',Best,'Worst',Worst,'FeasRate',FeasRate, ...
    'MeanViol',MeanViol,'MeanTime',MeanTime,'rank',rk,'avgRank',avgRank, ...
    'markSRDO',MK6,'pSRDO',P6,'a12SRDO',A12_6,'markNOLR',MK0,'pNOLR',P0,'a12NOLR',A12_0, ...
    'logFEs',logGrid,'medCurve',{medCurve},'allCurve',{allCurve}, ...
    'bestDesign',{bestDesign});
save(fullfile(cfg.data_dir,'analysis_eng.mat'),'A','-v7');
fprintf('[Statistics] Already written %s\n', xlsx);
end

% -------------------------------
function [MK,PM,AA] = pairtests(F, refName, algos)
[nP,nA,R]=size(F);
MK=strings(nP,nA); PM=nan(nP,nA); AA=nan(nP,nA);
base=find(strcmp(algos,refName),1);
for ip=1:nP
    ps=[]; pairs=[];
    for ai=1:nA
        if ai==base, continue; end
        x=squeeze(F(ip,base,:)); y=squeeze(F(ip,ai,:));
        ok=~isnan(x)&~isnan(y);
        AA(ip,ai)=A12(x(ok),y(ok));
        if sum(ok)>=6 && exist('signrank','file')==2
            p=signrank(x(ok),y(ok));
        else, p=NaN; end
        PM(ip,ai)=p; ps(end+1)=p; pairs(end+1)=ai; 
    end
    keep=~isnan(ps);
    padj=nan(size(ps));
    if any(keep), [~,~,pa]=holm(ps(keep)); padj(keep)=pa; end
    for k=1:numel(pairs)
        ai=pairs(k);
        x=squeeze(F(ip,base,:)); y=squeeze(F(ip,ai,:));
        if isnan(padj(k)), MK(ip,ai)='n/a';
        elseif padj(k)>=0.05, MK(ip,ai)='=';
        elseif nanmean(x)<nanmean(y), MK(ip,ai)='+';
        else, MK(ip,ai)='-'; end
    end
    MK(ip,base)='-';
end
end

function TW=wtbl(pnames,algos,MK,PM,AA,tag)
TW=table; TW.Problem=pnames(:);
for ai=1:numel(algos)
    TW.([algos{ai} '_mark'])=cellstr(MK(:,ai));
    TW.([algos{ai} '_p'])=PM(:,ai);
    TW.([algos{ai} '_A12'])=AA(:,ai);
end
TW.Properties.Description=sprintf('mark: + %s Significantly superior / = indifference / - %s Significantly inferior（Wilcoxon+Holm, alpha=0.05）',tag,tag);
end

function r=tr(x)
if exist('tiedrank','file')==2, r=tiedrank(x); else
    [~,ix]=sort(x); r=zeros(size(x)); r(ix)=1:numel(x);
    u=unique(x); for i=1:numel(u), r(x==u(i))=mean(r(x==u(i))); end
end
end
function a=A12(x,y)
x=x(:); y=y(:);
if isempty(x)||isempty(y), a=NaN; return; end
m=numel(x); n=numel(y);
a=(sum(x>y)+0.5*sum(x==y))/(m*n);
end
function [h,padj,praw]=holm(p)
m=numel(p); praw=p(:); [ps,ix]=sort(praw); adj=ps.*(m-(1:m).'); adj=min(adj,1);
for i=2:m, adj(i)=max(adj(i),adj(i-1)); end
padj=nan(m,1); padj(ix)=adj; h=padj<0.05;
end
