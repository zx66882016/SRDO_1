function fig_ablations(cfg)


d=fullfile(cfg.data_dir,'perrun_abl');
if isempty(dir(fullfile(d,'*.mat'))), warning('No ablation results, skip ablation image'); return; end

getlogs=@(group,tag,f) load_abl(d,group,tag,f);

% --------
vc=[0.70 0.70 0.70; 0.34 0.71 0.91; 0.00 0.62 0.45; 0.90 0.60 0.00; 0 0 0];
lss={':','--','--','-.','-'}; lws=[1.2 1.2 1.2 1.2 2.0];
fig=figure('Position',[50 50 720 620]);
tiledlayout(fig,ceil(numel(cfg.repFunc)/3),3,'Padding','compact','TileSpacing','compact');
for k=1:numel(cfg.repFunc)
    f=cfg.repFunc(k); ax=nexttile; hold(ax,'on');
    for v=0:4
        [x,Y]=getlogs('oper',sprintf('Oper_V%d',v),f);
        if isempty(x), continue; end
        semilogy(ax,x,max(median(Y,1),1e-12),'Color',vc(v+1,:),'LineStyle',lss{v+1},...
            'LineWidth',lws(v+1),'DisplayName',sprintf('O%d',v));
    end
    title(ax,sprintf('F%d',f)); grid(ax,'on'); xlabel(ax,'Function evaluations');
    ylabel(ax,'Error f(x)-f*'); xlim(ax,[1 cfg.fastBudget]); pub_style(fig,ax);
    if k==1, legend(ax,'Location','northeast'); end
end
pub_export(fig,'abl_operator_conv',cfg); close(fig);

% 
if exist('boxplot','file')==2
    fig=figure('Position',[50 50 720 620]);
    tiledlayout(fig,ceil(numel(cfg.repFunc)/3),3,'Padding','compact','TileSpacing','compact');
    for k=1:numel(cfg.repFunc)
        f=cfg.repFunc(k); ax=nexttile; M=[];
        for v=0:4
            [~,Y]=getlogs('oper',sprintf('Oper_V%d',v),f);
            if isempty(Y), M=[M nan(cfg.runs,1)]; else, M=[M max(Y(:,end),1e-12)]; end 
        end
        boxplot(ax,max(M,1e-12),'Labels',arrayfun(@(v)sprintf('O%d',v),0:4,'UniformOutput',false),...
            'Colors',[0.1 0.1 0.1],'Width',0.55);
        set(ax,'YScale','log'); title(ax,sprintf('F%d',f)); grid(ax,'on'); pub_style(fig,ax);
    end
    pub_export(fig,'abl_operator_box',cfg); close(fig);
end

% -------
sc=[0 0 0; 0.00 0.45 0.70; 0.84 0.37 0.00; 0.80 0.47 0.65];
lss={'-','--','-.',':'}; lws=[2 1.3 1.3 1.3];
fig=figure('Position',[50 50 720 620]);
tiledlayout(fig,ceil(numel(cfg.repFunc)/3),3,'Padding','compact','TileSpacing','compact');
for k=1:numel(cfg.repFunc)
    f=cfg.repFunc(k); ax=nexttile; hold(ax,'on');
    for s=1:4
        [x,Y]=getlogs('param',sprintf('Param_S%d_T500',s),f);
        if isempty(x), continue; end
        semilogy(ax,x,max(median(Y,1),1e-12),'Color',sc(s,:),'LineStyle',lss{s},...
            'LineWidth',lws(s),'DisplayName',sprintf('S%d',s));
    end
    title(ax,sprintf('F%d',f)); grid(ax,'on'); xlabel(ax,'Function evaluations');
    ylabel(ax,'Error f(x)-f*'); xlim(ax,[1 cfg.fastBudget]); pub_style(fig,ax);
    if k==1, legend(ax,'Location','northeast'); end
end
pub_export(fig,'abl_param_conv',cfg); close(fig);

% -------
rc=[0 0 0; 0.00 0.45 0.70; 0.84 0.37 0.00; 0.00 0.62 0.45];
lss={'-','--','-.',':'}; lws=[2 1.3 1.3 1.3];
fig=figure('Position',[50 50 720 250]);
tiledlayout(fig,1,numel(cfg.replaceFunc),'Padding','compact','TileSpacing','compact');
for k=1:numel(cfg.replaceFunc)
    f=cfg.replaceFunc(k); ax=nexttile; hold(ax,'on');
    for rp=0:3
        [x,Y]=getlogs('repl',sprintf('Repl_R%d',rp),f);
        if isempty(x), continue; end
        semilogy(ax,x,max(median(Y,1),1e-12),'Color',rc(rp+1,:),'LineStyle',lss{rp+1},...
            'LineWidth',lws(rp+1),'DisplayName',sprintf('R%d',rp));
    end
    title(ax,sprintf('F%d',f)); grid(ax,'on'); xlabel(ax,'Function evaluations');
    ylabel(ax,'Error f(x)-f*'); xlim(ax,[1 cfg.fastBudget]); pub_style(fig,ax);
    if k==1, legend(ax,'Location','northeast'); end
end
pub_export(fig,'abl_replace_conv',cfg); close(fig);
end

function [x,Y]=load_abl(d,group,tag,f)
x=[]; Y=[];
files=dir(fullfile(d,'*.mat'));
for i=1:numel(files)
    S=load(fullfile(d,files(i).name)); q=S.R;
    if ~isfield(q,'group'), continue; end
    if strcmp(q.group,group)&&strcmp(q.tag,tag)&&q.func==f
        if isempty(x), x=q.logFEs; end
        Y(end+1,:)=q.logError; 
    end
end
end
