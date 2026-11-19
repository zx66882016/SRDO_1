function fig_main(cfg)
for s=1:numel(cfg.run_suites)
    sname=cfg.run_suites{s};
    si=find(strcmp({cfg.suites.name},sname),1);
    S=cfg.suites(si);
    L=load_main(cfg,sname,S.stdBudget);
    Ll=load_main(cfg,sname,cfg.fastBudget);
    if isempty(fieldnames(L)), continue; end
    algos=fieldnames(L);
    order=cfg.main_algos; if cfg.include_legacy, order=[order,cfg.legacy_algos]; end
    algos=order(ismember(order,algos));

    repf = cfg.repFunc; repf=repf(ismember(repf,S.funcs));
    if strcmp(sname,'CEC2022'), repf=[1 3 6 8 10 12]; repf=repf(ismember(repf,S.funcs)); end
    if numel(repf)>6, repf=repf(1:6); end

    % --------------------
    fig=figure('Position',[50 50 720 620]);
    tiledlayout(fig,ceil(numel(repf)/3),3,'Padding','compact','TileSpacing','compact');
    for k=1:numel(repf)
        f=repf(k); ax=nexttile; hold(ax,'on');
        for ai=1:numel(algos)
            D=L.(algos{ai}); fi=find(D.funcs==f,1); if isempty(fi), continue; end
            x=D.logFEs{fi}; y=max(median(D.logs{fi},1),1e-12);
            st=algo_style(algos{ai});
            semilogy(ax,x,y,'Color',st.color,'LineStyle',st.ls,'LineWidth',1.4,...
                'Marker',st.mk,'MarkerIndices',round(linspace(1,numel(x),7)),...
                'MarkerSize',4,'DisplayName',st.disp);
        end
        title(ax,sprintf('%s F%d',sname,f)); grid(ax,'on');
        xlabel(ax,'Function evaluations'); ylabel(ax,'Error f(x)-f*');
        xlim(ax,[1 S.stdBudget]); pub_style(fig,ax);
        if k==1, legend(ax,'Location','northeast'); end
    end
    pub_export(fig,sprintf('conv_%s_rep',sname),cfg); close(fig);

    % ---------- ----------
    nf=numel(S.funcs); nc=5; nr=ceil(nf/nc);
    fig=figure('Position',[20 20 220*nc 180*nr]);
    tl=tiledlayout(fig,nr,nc,'Padding','compact','TileSpacing','compact'); ax=[];
    for k=1:nf
        f=S.funcs(k); ax=nexttile; hold(ax,'on');
        for ai=1:numel(algos)
            D=L.(algos{ai}); fi=find(D.funcs==f,1); if isempty(fi), continue; end
            x=D.logFEs{fi}; y=max(median(D.logs{fi},1),1e-12);
            st=algo_style(algos{ai});
            semilogy(ax,x,y,'Color',st.color,'LineStyle',st.ls,'LineWidth',1.1,...
                'DisplayName',st.disp);
        end
        title(ax,sprintf('F%d',f)); grid(ax,'on'); xlim(ax,[1 S.stdBudget]); pub_style(fig,ax);
    end
    lg=legend(ax,'Orientation','horizontal'); lg.Layout.Tile='south';
    xlabel(tl,'Function evaluations'); ylabel(tl,'Error f(x)-f*');
    pub_export(fig,sprintf('conv_%s_all',sname),cfg); close(fig);

    % --------------------
    cps=S.checkpoints;
    Rk=scanRanks(L,algos,S.funcs,cps,cfg.runs);
    fig=figure('Position',[50 50 460 330]); ax=axes(fig); hold(ax,'on');
    for ai=1:numel(algos)
        st=algo_style(algos{ai});
        plot(ax,cps,Rk(:,ai),'Color',st.color,'LineStyle',st.ls,'LineWidth',1.8,...
            'Marker',st.mk,'MarkerSize',5,'DisplayName',st.disp);
    end
    set(ax,'XScale','log'); grid(ax,'on'); axis(ax,'tight');
    xlabel(ax,'Function evaluations'); ylabel(ax,'Average Friedman rank');
    title(ax,sname);
    leg=legend(ax,'Location','southoutside','Orientation','horizontal');
    leg.FontSize=8; leg.Box='off';
    pub_style(fig,ax);
    pub_export(fig,sprintf('rank_scan_%s',sname),cfg); close(fig);

    % ------------------
    rl=scanRanks(Ll,algos,S.funcs,cfg.fastBudget,cfg.runs);
    rs=scanRanks(L, algos,S.funcs,S.stdBudget,cfg.runs);
    fig=figure('Position',[50 50 460 300]); ax=axes(fig);
    b=bar(ax,[rl; rs]');
    b(1).FaceColor=[0.85 0.85 0.85]; b(2).FaceColor=[0.35 0.35 0.35];
    set(ax,'XTick',1:numel(algos),'XTickLabel',cellfun(@(a)algo_style(a).disp,algos,'UniformOutput',false));
    ylabel(ax,'Average Friedman rank'); set(ax,'YGrid','on','XGrid','off');
    legend(ax,{sprintf('Limited (%g FEs)',cfg.fastBudget),sprintf('Standard (%g FEs)',S.stdBudget)},...
           'Location','northwest');
    pub_style(fig,ax); pub_export(fig,sprintf('rank_bar_%s',sname),cfg); close(fig);
end
end

% ----------------aveage rank---------------
function Rk=scanRanks(L,algos,funcs,cps,runs)
cps=cps(:).';
Rk=nan(numel(cps),numel(algos));
for ci=1:numel(cps)
    per=nan(numel(funcs)*runs,numel(algos)); row=0;
    for k=1:numel(funcs)
        D0=cell(1,numel(algos)); ok=true;
        for ai=1:numel(algos)
            if ~isfield(L,algos{ai}), ok=false; break; end
            D=L.(algos{ai}); fi=find(D.funcs==funcs(k),1);
            if isempty(fi), ok=false; break; end
            ix=find(D.cpFEs{fi}>=cps(ci),1);
            if isempty(ix), ok=false; break; end
            D0{ai}=D.cpErr{fi}(:,ix);
        end
        if ~ok, continue; end
        for r=1:runs
            row=row+1; v=cellfun(@(z)z(r),D0);
            per(row,:)=tiedrank_local(v);
        end
    end
    Rk(ci,:)=mean(per(1:row,:),1,'omitnan');
end
end

% ----------------data---------------
function L=load_main(cfg,sname,budget)
L=struct();
d=fullfile(cfg.data_dir,sprintf('perrun_main_B%d',budget));
files=dir(fullfile(d,[sname '_*.mat']));
for i=1:numel(files)
    S=load(fullfile(d,files(i).name)); q=S.R;
    if ~isfield(L,q.algo), L.(q.algo)=struct('funcs',[],'cpFEs',{{}},'cpErr',{{}},'logFEs',{{}},'logs',{{}}); end
    D=L.(q.algo); fi=find(D.funcs==q.func,1);
    if isempty(fi)
        fi=numel(D.funcs)+1; D.funcs(fi)=q.func;
        D.cpFEs{fi}=q.cpFEs; D.logFEs{fi}=q.logFEs;
        D.cpErr{fi}=zeros(0,numel(q.cpError));
        D.logs{fi}=zeros(0,numel(q.logError));
    end
    D.cpErr{fi}(end+1,:)=q.cpError(:).';
    D.logs{fi}(end+1,:)=q.logError(:).';
    L.(q.algo)=D;
end
end
function r=tiedrank_local(v)
[~,ix]=sort(v); r=zeros(size(v)); r(ix)=1:numel(v);
u=unique(v);
for i=1:numel(u), r(v==u(i))=mean(r(v==u(i))); end
end
