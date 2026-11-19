function pub_style(fig, ax)
if nargin<2 || isempty(ax), ax=gca(fig); end
set(0,'defaultAxesFontName','Arial');
set(fig,'Color','w','Units','inches');
set(ax,'FontName','Arial','FontSize',9,'TickDir','in','Box','on',...
       'XMinorTick','on','YMinorTick','on','LineWidth',0.8,...
       'XColor',[0.1 0.1 0.1],'YColor',[0.1 0.1 0.1],...
       'GridColor',[0.85 0.85 0.85],'GridAlpha',0.6,'MinorGridAlpha',0.25);
set(ax.XLabel,'FontName','Arial','FontSize',10);
set(ax.YLabel,'FontName','Arial','FontSize',10);
set(ax.Title,'FontName','Arial','FontSize',10,'FontWeight','normal');
leg=findobj(fig,'Type','legend');
for k=1:numel(leg)
    set(leg(k),'FontName','Arial','FontSize',8,'Box','off','EdgeColor','none');
end
end

function c = wong_colors()
c = [0.00 0.00 0.00;   
     0.90 0.60 0.00;   
     0.34 0.71 0.91;   
     0.00 0.62 0.45;   
     0.94 0.89 0.26;   
     0.00 0.45 0.70;   
     0.84 0.37 0.00;   
     0.80 0.47 0.65];  
end
