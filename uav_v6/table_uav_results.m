function table_uav_results(cfg)


if nargin<1||isempty(cfg), cfg=uav_config(); end
S=load(fullfile(cfg.data_dir,'uav_results.mat'));
data24=S.data24; data48=S.data48; astar=S.astar;
nS=numel(cfg.scenarios); nA=size(cfg.algos,1);
keys=cfg.algos(:,1); disps=cfg.algos(:,2);
out=fullfile(cfg.table_dir,'uav_results.xlsx');
if exist(out,'file'), delete(out); end


summarizeBlock(data24, nS, nA, cfg, keys, disps, out, 'Main_24', true);
block48=cell(1,nA,cfg.runs);
for ai=1:nA, for r=1:cfg.runs, block48{1,ai,r}=data48{ai,r}; end, end
summarizeBlock(block48, 1, nA, cfg, keys, disps, out, 'HighDim_48', false);

%% ---- A* ----
rows=[];
for s=1:nS
    a=astar{s}; sc=uav_scenarios(cfg.scenarios{s});
    if isempty(a)||~isfield(a,'found')||~a.found
        rows=[rows; {sc.label,'A* raw','not found',nan,nan,nan,nan,nan,a.searchTime}]; 
        rows=[rows; {sc.label,'A* + spline','not found',nan,nan,nan,nan,nan,nan}]; 
        continue;
    end
    rw=a.raw; sp=a.spline;
    rows=[rows; {sc.label,'A* raw','feasible',rw.cost,rw.len,rw.curv, ...
        double(rw.feasible),rw.min_clearance,a.searchTime}]; 
    rows=[rows; {sc.label,'A* + spline',feasWord(sp.feasible),sp.cost,sp.len,sp.curv, ...
        double(sp.feasible),sp.min_clearance,a.splineTime}]; 
end
T=cell2table(rows,'VariableNames',{'Scenario','Variant','Status','Cost', ...
    'PathLength','Curvature','Feasible','MinClearance','Time_s'});
writetable(T,out,'Sheet','Astar');

%% ---- Wilcoxon----
wrows=[];
aiSR=find(strcmp(keys,'SRDO'));
for s=1:nS
    base=arrayfun(@(r)data24{s,aiSR,r}.best,1:cfg.runs);
    pvals=nan(1,nA); marks=repmat({'='},1,nA);
    for ai=1:nA
        if ai==aiSR, continue; end
        comp=arrayfun(@(r)data24{s,ai,r}.best,1:cfg.runs);
        pvals(ai)=signrank(base,comp);
    end
    pv=pvals; pv(isnan(pv))=1;
    [~,padj,~]=fdr_holm(pv,0.05);
    for ai=1:nA
        if ai==aiSR, continue; end
        comp=arrayfun(@(r)data24{s,ai,r}.best,1:cfg.runs);
        if padj(ai)>=0.05, marks{ai}='=';
        elseif median(base)<median(comp), marks{ai}='+';   
        else, marks{ai}='-';
        end
        wrows=[wrows; {uav_scenarios(cfg.scenarios{s}).label,disps{ai}, ...
            marks{ai},pvals(ai),padj(ai)}]; 
    end
end
%  S3
base=arrayfun(@(r)data48{aiSR,r}.best,1:cfg.runs);
pvals=nan(1,nA);
for ai=1:nA
    if ai==aiSR, continue; end
    comp=arrayfun(@(r)data48{ai,r}.best,1:cfg.runs);
    pvals(ai)=signrank(base,comp);
end
pv=pvals; pv(isnan(pv))=1; [~,padj,~]=fdr_holm(pv,0.05);
for ai=1:nA
    if ai==aiSR, continue; end
    comp=arrayfun(@(r)data48{ai,r}.best,1:cfg.runs);
    if padj(ai)>=0.05, mk='='; elseif median(base)<median(comp), mk='+'; else, mk='-'; end
    wrows=[wrows; {'S3 dense D=48',disps{ai},mk,pvals(ai),padj(ai)}]; 
end
T=cell2table(wrows,'VariableNames',{'Scenario','Comparator','Mark','p_raw','p_Holm'});
writetable(T,out,'Sheet','Wilcoxon');
fprintf('[Table] Exported %s\n',out);
end

% ---------------- ---------------
function summarizeBlock(block,nS,nA,cfg,keys,disps,out,sheetName,multi)
rows=[];
for s=1:nS
    means=nan(1,nA);
    for ai=1:nA
        vals=arrayfun(@(r)block{s,ai,r}.best,1:cfg.runs);
        feas=arrayfun(@(r)block{s,ai,r}.metrics.feasible,1:cfg.runs);
        lens=arrayfun(@(r)block{s,ai,r}.metrics.len,1:cfg.runs);
        tims=arrayfun(@(r)block{s,ai,r}.time,1:cfg.runs);
        means(ai)=mean(vals);
        if multi, lab=uav_scenarios(cfg.scenarios{s}).label; else, lab='S3 dense (16 pts)'; end
        rows=[rows; {lab,disps{ai},min(vals),max(vals),mean(vals),std(vals), ...
            median(vals),mean(feas),mean(lens),mean(tims),nan}]; 
    end
    % 
    idx=find(strcmp(rows(:,1),lab));
    [~,ord]=sort(means); rk=zeros(1,nA); rk(ord)=1:nA;
    for j=1:numel(idx), rows{idx(j),11}=rk(j); end
end
T=cell2table(rows,'VariableNames',{'Scenario','Algorithm','Best','Worst', ...
    'Mean','Std','Median','FeasRate','MeanPathLen','MeanTime_s','Rank'});
writetable(T,out,'Sheet',sheetName);
end

function w=feasWord(f)
if f, w='feasible'; else, w='infeasible'; end
end

function [reject,padj,zcrit]=fdr_holm(p,alpha)
% Holm 
m=numel(p); [ps,ord]=sort(p(:)');
padj_sorted=nan(1,m);
running=0;
for i=1:m
    running=max(running,(m-i+1)*ps(i));
    padj_sorted(i)=min(1,running);
end
padj=nan(1,m); padj(ord)=padj_sorted;
reject=padj<alpha; zcrit=nan;
end
