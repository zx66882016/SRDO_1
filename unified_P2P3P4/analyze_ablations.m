function AA = analyze_ablations(cfg)
outdir=fullfile(cfg.data_dir,'perrun_abl');
files=dir(fullfile(outdir,'*.mat'));
if isempty(files), warning('Ablation results not found'); AA=struct(); return; end
G=struct();
for i=1:numel(files)
    S=load(fullfile(outdir,files(i).name)); q=S.R;
    if ~isfield(q,'group')||~ismember(q.group,{'param','oper','repl'}), continue; end
    key=sprintf('%s__%s__F%d',q.group,q.tag,q.func);
    if ~isfield(G,mk(key)), G.(mk(key))=[]; end
    G.(mk(key))(end+1)=q.finalError;
end
xlsx=fullfile(cfg.table_dir,'ablation_results.xlsx');
if exist(xlsx,'file'), delete(xlsx); end
AA=struct();

% ---- 1) Parameter Ablation ----
for bud={'T500','T200'}
    tags=arrayfun(@(s) sprintf('Param_S%d_%s',s,bud{1}),1:4,'UniformOutput',false);
    funcs=cfg.repFunc;
    M=nan(numel(funcs),4); SD=M; MN=M;
    P=nan(numel(funcs),3); MK=strings(numel(funcs),3);
    for fi=1:numel(funcs)
        vals=[];
        for s=1:4
            key=sprintf('param__%s__F%d',tags{s},funcs(fi));
            if isfield(G,mk(key)), v=G.(mk(key)); else, v=[]; end
            M(fi,s)=mean(v); SD(fi,s)=std(v); MN(fi,s)=min(v); vals(s,1:numel(v))=v;
        end
        if size(vals,2)>=5
            ps=arrayfun(@(s) safe_signrank(vals(1,:),vals(s,:)),2:4);
            [~,~,pa]=holm(ps); P(fi,:)=pa;
            for s=1:3, MK(fi,s)=mkmark(pa(s),mean(vals(1,:)),mean(vals(s+1,:))); end
        end
    end
    T=table;
    T.Func=arrayfun(@(f)sprintf('F%d',f),funcs,'UniformOutput',false).';
    for s=1:4
        T.(sprintf('S%d_mean',s))=M(:,s); T.(sprintf('S%d_std',s))=SD(:,s); T.(sprintf('S%d_min',s))=MN(:,s);
    end
    T.p_S2vsS1=P(:,1); T.mark_S2=cellstr(MK(:,1));
    T.p_S3vsS1=P(:,2); T.mark_S3=cellstr(MK(:,2));
    T.p_S4vsS1=P(:,3); T.mark_S4=cellstr(MK(:,3));
    writetable(T,xlsx,'Sheet',['param_' bud{1}],'WriteMode','overwritesheet');
    AA.(['param_' bud{1}])=struct('funcs',funcs,'M',M,'SD',SD,'MN',MN,'P',P,'MK',MK);
end

% ---- 2)V0..V4 ----
funcs=cfg.repFunc;
M=nan(numel(funcs),5); SD=M; MN=M; P=nan(numel(funcs),4); MK=strings(numel(funcs),4);
for fi=1:numel(funcs)
    vals=[];
    for v=0:4
        key=sprintf('oper__Oper_V%d__F%d',v,funcs(fi));
        if isfield(G,mk(key)), vv=G.(mk(key)); else, vv=[]; end
        M(fi,v+1)=mean(vv); SD(fi,v+1)=std(vv); MN(fi,v+1)=min(vv); vals(v+1,1:numel(vv))=vv;
    end
    if size(vals,2)>=5
        ps=arrayfun(@(v) safe_signrank(vals(5,:),vals(v,:)),1:4); % V4 vs V0..V3
        [~,~,pa]=holm(ps); P(fi,:)=pa;
        for v=1:4, MK(fi,v)=mkmark(pa(v),mean(vals(5,:)),mean(vals(v,:))); end
    end
end
T=table; T.Func=arrayfun(@(f)sprintf('F%d',f),funcs,'UniformOutput',false).';
for v=0:4
    T.(sprintf('V%d_mean',v))=M(:,v+1); T.(sprintf('V%d_std',v))=SD(:,v+1); T.(sprintf('V%d_min',v))=MN(:,v+1);
end
for v=1:4
    T.(sprintf('p_V4vsV%d',v-1))=P(:,v); T.(sprintf('mark_V%d',v-1))=cellstr(MK(:,v));
end
writetable(T,xlsx,'Sheet','operator','WriteMode','overwritesheet');
AA.operator=struct('funcs',funcs,'M',M,'SD',SD,'MN',MN,'P',P,'MK',MK);

% ---- 3) R0..R3 ----
funcs=cfg.replaceFunc;
M=nan(numel(funcs),4); SD=M; MN=M; P=nan(numel(funcs),3); MK=strings(numel(funcs),3);
for fi=1:numel(funcs)
    vals=[];
    for rp=0:3
        key=sprintf('repl__Repl_R%d__F%d',rp,funcs(fi));
        if isfield(G,mk(key)), vv=G.(mk(key)); else, vv=[]; end
        M(fi,rp+1)=mean(vv); SD(fi,rp+1)=std(vv); MN(fi,rp+1)=min(vv); vals(rp+1,1:numel(vv))=vv;
    end
    if size(vals,2)>=5
        ps=arrayfun(@(rp) safe_signrank(vals(1,:),vals(rp,:)),2:4);
        [~,~,pa]=holm(ps); P(fi,:)=pa;
        for rp=1:3, MK(fi,rp)=mkmark(pa(rp),mean(vals(1,:)),mean(vals(rp+1,:))); end
    end
end
T=table; T.Func=arrayfun(@(f)sprintf('F%d',f),funcs,'UniformOutput',false).';
for rp=0:3
    T.(sprintf('R%d_mean',rp))=M(:,rp+1); T.(sprintf('R%d_std',rp))=SD(:,rp+1); T.(sprintf('R%d_min',rp))=MN(:,rp+1);
end
T.p_R1=P(:,1); T.mark_R1=cellstr(MK(:,1));
T.p_R2=P(:,2); T.mark_R2=cellstr(MK(:,2));
T.p_R3=P(:,3); T.mark_R3=cellstr(MK(:,3));
writetable(T,xlsx,'Sheet','replacement','WriteMode','overwritesheet');
AA.replacement=struct('funcs',funcs,'M',M,'SD',SD,'MN',MN,'P',P,'MK',MK);

save(fullfile(cfg.data_dir,'analysis_abl.mat'),'AA','-v7');
fprintf('[Ablation statistics] Already written %s\n',xlsx);
fprintf('Marking Instructions：+ Complete alignment is significantly better；- The variant is significantly superior；= No significant difference（Holm corrected p<0.05）\n');
end
function s=mk(k), s=['k' char(regexprep(k,'[^a-zA-Z0-9]','_'))]; end
function p=safe_signrank(a,b)
if exist('signrank','file')==2 && numel(a)==numel(b), p=signrank(a,b); else, p=NaN; end
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