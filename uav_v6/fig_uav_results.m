function fig_uav_results(cfg)

if nargin<1||isempty(cfg), cfg=uav_config(); end
outFile=fullfile(cfg.data_dir,'uav_results.mat');
S=load(outFile);
data24=S.data24; data48=S.data48; astar=S.astar;
nS=numel(cfg.scenarios); nA=size(cfg.algos,1);
keys=cfg.algos(:,1);

%% 
fig=figure('Position',[40 40 900 720]);
tiledlayout(fig,2,2,'Padding','compact','TileSpacing','compact');
for s=1:nS
    ax=nexttile; hold(ax,'on');
    B=cfg.fes24; gridFes=unique(round(logspace(0,log10(B),200)));
    for ai=1:nA
        Y=nan(cfg.runs,numel(gridFes));
        for r=1:cfg.runs
            R=data24{s,ai,r};
            Y(r,:)=interp1(R.logFEs,R.logBest,gridFes,'previous','extrap');
        end
        st=uav_style(keys{ai});
        plot(ax,gridFes,max(median(Y,1),1e-8),'Color',st.color, ...
            'LineStyle',st.ls,'LineWidth',st.lw,'Marker',st.mk, ...
            'MarkerSize',4,'MarkerIndices',round(linspace(1,numel(gridFes),7)), ...
            'DisplayName',st.disp);
    end
    set(ax,'YScale','log','XScale','log'); grid(ax,'on');
    xlabel(ax,'Function evaluations'); ylabel(ax,'Cost');
    title(ax,uav_scenarios(cfg.scenarios{s}).label); xlim(ax,[cfg.N cfg.fes24]);
    if s==1, legend(ax,'Location','northeast','FontSize',7,'NumColumns',2); end
    pub_style(fig,ax);
end
pub_export(fig,'uav_conv_24',cfg); close(fig);

%% 
if exist('boxplot','file')==2
    fig=figure('Position',[40 40 900 720]);
    tiledlayout(fig,2,2,'Padding','compact','TileSpacing','compact');
    for s=1:nS
        ax=nexttile;
        M=nan(cfg.runs,nA);
        for ai=1:nA
            for r=1:cfg.runs, M(r,ai)=max(data24{s,ai,r}.best,1e-8); end
        end
        labs=cellfun(@(k)shortlab(k),keys,'UniformOutput',false);
        boxplot(ax,M,'Labels',labs,'Colors',[0.1 0.1 0.1],'Width',0.6,'LabelOrientation','inline');
        set(ax,'YScale','log'); grid(ax,'on');
        title(ax,uav_scenarios(cfg.scenarios{s}).label);
        ylabel(ax,'Final cost');
        pub_style(fig,ax);
    end
    pub_export(fig,'uav_box_24',cfg); close(fig);
end

%% 
fig=figure('Position',[40 40 1000 800]);
tiledlayout(fig,2,2,'Padding','compact','TileSpacing','compact');
for s=1:nS
    ax=nexttile; hold(ax,'on');
    sc=uav_scenarios(cfg.scenarios{s});
    drawObstacles(ax,sc.obs);
    %
    aiSR=find(strcmp(keys,'SRDO'));
    best=inf; bestPath=[];
    for r=1:cfg.runs
        R=data24{s,aiSR,r};
        if R.metrics.feasible && R.best<best, best=R.best; bestPath=R.metrics.path; end
    end
    if isempty(bestPath)  
        for r=1:cfg.runs
            R=data24{s,aiSR,r};
            if R.best<best, best=R.best; bestPath=R.metrics.path; end
        end
    end
    plot3(ax,bestPath(:,1),bestPath(:,2),bestPath(:,3),'-', ...
        'Color',[0.80 0.15 0.15],'LineWidth',2.2,'DisplayName','SRDO');
  
    if ~isempty(astar{s}) && astar{s}.found && astar{s}.spline.feasible
        pa=astar{s}.spline.path;
        plot3(ax,pa(:,1),pa(:,2),pa(:,3),'--','Color',[0.0 0.35 0.65], ...
            'LineWidth',1.8,'DisplayName','A* + spline');
    end
    plot3(ax,sc.start(1),sc.start(2),sc.start(3),'go','MarkerFaceColor','g','MarkerSize',8,'DisplayName','Start');
    plot3(ax,sc.goal(1),sc.goal(2),sc.goal(3),'rs','MarkerFaceColor','r','MarkerSize',8,'DisplayName','Goal');
    ax.XLim=[0 100];ax.YLim=[0 100];ax.ZLim=[0 100];
    xlabel(ax,'X (m)');ylabel(ax,'Y (m)');zlabel(ax,'Z (m)');
    title(ax,sc.label); grid(ax,'on'); view(ax,38,28); axis(ax,'equal');
    if s==1, legend(ax,'Location','northeast','FontSize',8); end
    pub_style(fig,ax);
end
pub_export(fig,'uav_paths_3d',cfg); close(fig);


fig=figure('Position',[40 40 1000 420]);
tiledlayout(fig,1,2,'Padding','compact','TileSpacing','compact');
ax=nexttile; hold(ax,'on');
B=cfg.fes48; gridFes=unique(round(logspace(0,log10(B),200)));
for ai=1:nA
    Y=nan(cfg.runs,numel(gridFes));
    for r=1:cfg.runs
        R=data48{ai,r};
        Y(r,:)=interp1(R.logFEs,R.logBest,gridFes,'previous','extrap');
    end
    st=uav_style(keys{ai});
    plot(ax,gridFes,max(median(Y,1),1e-8),'Color',st.color,'LineStyle',st.ls, ...
        'LineWidth',st.lw,'Marker',st.mk,'MarkerSize',4, ...
        'MarkerIndices',round(linspace(1,numel(gridFes),7)),'DisplayName',st.disp);
end
set(ax,'YScale','log','XScale','log'); grid(ax,'on');
xlabel(ax,'Function evaluations'); ylabel(ax,'Cost');
title(ax,'S3 dense, 16 control points (D=48)');
legend(ax,'Location','northeast','FontSize',7,'NumColumns',2); xlim(ax,[cfg.N cfg.fes48]);
pub_style(fig,ax);

if exist('boxplot','file')==2
    ax=nexttile;
    M=nan(cfg.runs,nA);
    for ai=1:nA
        for r=1:cfg.runs, M(r,ai)=max(data48{ai,r}.best,1e-8); end
    end
    labs=cellfun(@(k)shortlab(k),keys,'UniformOutput',false);
    boxplot(ax,M,'Labels',labs,'Colors',[0.1 0.1 0.1],'Width',0.6,'LabelOrientation','inline');
    set(ax,'YScale','log'); grid(ax,'on');
    ylabel(ax,'Final cost'); title(ax,sprintf('S3 dense, D=48, %d runs',cfg.runs));
    pub_style(fig,ax);
end
pub_export(fig,'uav_highdim_48',cfg); close(fig);


ranks24=nan(nS,nA); ranks48=nan(1,nA);
for s=1:nS
    perRun=nan(cfg.runs,nA);
    for r=1:cfg.runs
        vals=arrayfun(@(ai)data24{s,ai,r}.best,1:nA);
        perRun(r,:)=tiedrank(vals);
    end
    ranks24(s,:)=mean(perRun,1);
end
for r=1:cfg.runs
    vals=arrayfun(@(ai)data48{ai,r}.best,1:nA);
    ranks48(r,:)=tiedrank(vals);
end
ranks48=mean(ranks48,1);
meanRank=mean(ranks24,1);
[~,ord]=sort(meanRank);
fig=figure('Position',[40 40 820 420]);
b=bar(meanRank(ord),'FaceColor',[0.45 0.65 0.85],'EdgeColor',[0.1 0.1 0.1]);
hold on;
plot(ranks48(ord),'d','Color',[0.80 0.15 0.15],'MarkerFaceColor',[0.80 0.15 0.15], ...
    'MarkerSize',8,'LineStyle','none','DisplayName','D=48 (S3)');
set(gca,'XTick',1:nA,'XTickLabel',cellfun(@(k)uav_style(k).disp,keys(ord),'UniformOutput',false));
xtickangle(25); ylabel('Average rank (lower is better)'); grid on;
legend({'D=24, mean over 4 scenarios','D=48, S3 dense'},'Location','northwest');
title('Average rank across UAV scenarios');
pub_export(fig,'uav_rank_bar',cfg); close(fig);
end


function st=uav_style(key)
switch key
    case 'SRDO',      st=algo_style('SRDO');       st.lw=2.0;
    case 'SRDO_NOLR', st=algo_style('SRDO_NOLR');  st.ls='--'; st.lw=1.4;
    otherwise,        st=algo_style(key);          st.lw=1.4;
end
end

function drawObstacles(ax,obs)
[xs,ys,zs]=sphere(14);
for k=1:size(obs,1)
    r=obs(k,4); c=obs(k,1:3);
    surf(ax,c(1)+r*xs,c(2)+r*ys,c(3)+r*zs, ...
        'FaceColor',[0.75 0.78 0.82],'EdgeColor','none', ...
        'FaceAlpha',0.35,'HandleVisibility','off');
end
end

function s=shortlab(key)

M=containers.Map({'SRDO','SRDO_NOLR','LSHADE','IPOP_CMAES','SFOA','GLO','PSO','GA','SAR','DBO','MSO','HBO','BNCO'}, ...
                 {'SRDO','w/o LR','L-SHADE','IPOP','SFOA','GLO','PSO','GA','SAR','DBO','MSO','HBO','BNCO'});
if isKey(M,key), s=M(key); else, s=key; end
end
