function fig_engineering(cfg)
if nargin < 1, cfg = eng_config(); end
S = load(fullfile(cfg.data_dir,'analysis_eng.mat')); A = S.A;
algos = A.algos; nP = numel(A.pids);
keyAlgos = cfg.key_algos;
keyAlgos = keyAlgos(ismember(keyAlgos,algos));

% --------------------------------
nameMap = containers.Map({ ...
    'SRDO','SRDO_NOLR','SRDO_HP','LSHADE','IPOP_CMAES', ...
    'PSO','GA','SAR','DBO','MSO','SFOA','GLO','HBO','BNCO'}, { ...
    'SRDO','SRDO (w/o LR)','SRDO-HP','L-SHADE','IPOP-CMA-ES', ...
    'PSO','GA','SAR','DBO','MSO','SFOA','GLO','HBO','BNCO'});
dn = @(a) nameMap(a);

% --------------------------------
probMap = containers.Map({ ...
    'welded_beam','pressure_vessel','speed_reducer', ...
    'three_bar_truss','rc_beam','tubular_column'}, { ...
    'Welded Beam','Pressure Vessel','Speed Reducer', ...
    'Three-bar Truss','Reinforced Concrete Beam','Tubular Column'});
pn = @(tag) probMap(tag);

labs   = cellfun(dn, keyAlgos, 'UniformOutput', false);
allLabs = cellfun(dn, algos,   'UniformOutput', false);

% --------------------------------
fig = figure('Position',[40 40 900 560],'MenuBar','none','ToolBar','none');
tiledlayout(fig,2,3,'Padding','compact','TileSpacing','compact');
for ip = 1:nP
    ax = nexttile; hold(ax,'on');
    for ai = 1:numel(keyAlgos)
        k = find(strcmp(algos,keyAlgos{ai}),1);
        if isempty(k), continue; end
        y = A.medCurve{ip,k}; x = A.logFEs;
        st = algo_style(keyAlgos{ai});
        n = numel(y); mi = unique(round(linspace(1,n,7)));
        semilogy(ax,x,y,'Color',st.color,'LineStyle',st.ls,'LineWidth',1.4, ...
            'Marker',st.mk,'MarkerIndices',mi,'MarkerSize',4.5, ...
            'DisplayName',dn(keyAlgos{ai}));
    end
    grid(ax,'on'); xlim(ax,[1 cfg.maxFEs]);
    title(ax,pn(A.pnames{ip}),'FontSize',11,'FontWeight','bold');
    xlabel(ax,'Function evaluations','FontSize',10);
    ylabel(ax,'Best feasible objective','FontSize',10);
    set(ax,'FontSize',9);
    if ip==1, legend(ax,'Location','northeast','FontSize',8); end
end
pub_export(fig,'eng_conv_key',cfg); close(fig);

% --------------------------------
fig = figure('Position',[40 40 900 560],'MenuBar','none','ToolBar','none');
tiledlayout(fig,2,3,'Padding','compact','TileSpacing','compact');
for ip = 1:nP
    ax = nexttile;
    M = nan(cfg.runs,numel(keyAlgos));
    for ai=1:numel(keyAlgos)
        k=find(strcmp(algos,keyAlgos{ai}),1);
        if ~isempty(k), M(:,ai)=squeeze(A.F(ip,k,:)); end
    end
    boxplot(ax,M,'Labels',labs,'Colors',[0.1 0.1 0.1],'Width',0.6);
    hb=findobj(ax,'Tag','Box');
    for i=1:numel(hb)
        idx=numel(hb)-i+1;
        if idx>=1 && idx<=numel(keyAlgos)
            st=algo_style(keyAlgos{idx});
            patch(get(hb(i),'XData'),get(hb(i),'YData'),st.color, ...
                'FaceAlpha',0.25,'EdgeColor',st.color,'Tag','Box');
        end
    end
    grid(ax,'on');
    title(ax,pn(A.pnames{ip}),'FontSize',11,'FontWeight','bold');
    ylabel(ax,'Final feasible objective','FontSize',10);
    set(ax,'FontSize',8);
    ax.XTickLabelRotation=30;
end
pub_export(fig,'eng_box_key',cfg); close(fig);

% ----------------+++++++++++----------------
[rkSort,ix]=sort(A.avgRank,'descend');
fig=figure('Position',[40 40 720 520],'MenuBar','none','ToolBar','none');
ax=axes(fig); hold(ax,'on');
cols=zeros(numel(ix),3);
for i=1:numel(ix), cols(i,:)=algo_style(algos{ix(i)}).color; end
b=barh(ax,rkSort,0.66,'FaceColor','flat'); b.CData=cols; b.FaceAlpha=.85;
for i=1:numel(ix)
    text(ax,rkSort(i)+0.08,i,sprintf('%.2f',rkSort(i)), ...
        'FontName','Arial','FontSize',9,'VerticalAlignment','middle');
end
set(ax,'YTick',1:numel(ix),'YTickLabel',cellfun(dn,algos(ix),'UniformOutput',false));
set(ax,'FontSize',10);
xlabel(ax,'Average Friedman rank over 6 problems (lower is better)','FontSize',11);
xlim(ax,[0 max(rkSort)+1.2]); grid(ax,'on');
pub_export(fig,'eng_rank_bar',cfg); close(fig);

% -------------------------------
for ip=1:nP
    fig=figure('Position',[40 40 680 480],'MenuBar','none','ToolBar','none');
    ax=axes(fig); hold(ax,'on');
    for ai=1:numel(algos)
        y=A.medCurve{ip,ai}; st=algo_style(algos{ai});
        n=numel(y); mi=unique(round(linspace(1,n,6)));
        semilogy(ax,A.logFEs,y,'Color',st.color,'LineStyle',st.ls,'LineWidth',1.1, ...
            'Marker',st.mk,'MarkerIndices',mi,'MarkerSize',4, ...
            'DisplayName',dn(algos{ai}));
    end
    grid(ax,'on'); xlim(ax,[1 cfg.maxFEs]);
    title(ax,pn(A.pnames{ip}),'FontSize',12,'FontWeight','bold');
    xlabel(ax,'Function evaluations','FontSize',10);
    ylabel(ax,'Best feasible objective','FontSize',10);
    set(ax,'FontSize',9);
    legend(ax,'Location','eastoutside','FontSize',7);
    pub_export(fig,sprintf('eng_conv_all_%s',A.pnames{ip}),cfg); close(fig);
end

% --------------------------------
fig=figure('Position',[40 40 1100 650],'MenuBar','none','ToolBar','none');
tiledlayout(fig,2,3,'Padding','compact','TileSpacing','compact');
for ip=1:nP
    ax=nexttile;
    M=squeeze(A.F(ip,:,:)).'; 
    boxplot(ax,M,'Labels',allLabs,'Colors',[0.1 0.1 0.1],'Width',0.6);
    hb=findobj(ax,'Tag','Box');
    for i=1:numel(hb)
        idx=numel(hb)-i+1;
        if idx>=1 && idx<=numel(algos)
            st=algo_style(algos{idx});
            patch(get(hb(i),'XData'),get(hb(i),'YData'),st.color, ...
                'FaceAlpha',0.2,'EdgeColor',st.color,'Tag','Box');
        end
    end
    grid(ax,'on');
    title(ax,pn(A.pnames{ip}),'FontSize',11,'FontWeight','bold');
    ylabel(ax,'Final feasible objective','FontSize',10);
    set(ax,'FontSize',7);
    ax.XTickLabelRotation=45;
end
pub_export(fig,'eng_box_all',cfg); close(fig);
fprintf('[Figure] The export of all engineering problem charts has been completed\n');
end
