function pub_export(fig, name, cfg)

out=cfg.fig_dir;
if ~exist(out,'dir'), mkdir(out); end
base=fullfile(out,name);
set(fig,'InvertHardcopy','off');
try, exportgraphics(fig,[base '.pdf'],'ContentType','vector'); catch ME, warning(ME.message); end
try, print(fig,[base '.svg'],'-dsvg','-vector'); catch ME, warning(ME.message); end
try, exportgraphics(fig,[base '.png'],'Resolution',600); catch ME, warning(ME.message); end
try, savefig(fig,[base '.fig']); catch ME, warning(ME.message); end
fprintf('[Figure] %s (.pdf/.svg/.png@600dpi/.fig)\n',name);
end
