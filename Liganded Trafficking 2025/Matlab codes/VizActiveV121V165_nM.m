function activeFig = VizActiveV121V165_nM(Routs,vals121,vals165,yaxistitle,figname,ngml,visibility,ligandmap)

v165ngmltonM = 1./44.;  % inverse of MW in kDa
v121ngmltonM = 1./28.;  % inverse of MW in kDa
pgf1ngmltonM = 1./29.7; % inverse of MW in kDa
pgf2ngmltonM = 1./34.6; % inverse of MW in kDa

minngml = min(ngml(2:end))*v121ngmltonM; % adjust to range encompassing all (hence V121 low)
maxngml = max(ngml(2:end))*v165ngmltonM; % adjust to range encompassing all (hence V165 high)

a121 = []; b121 = []; c121 = [];
a165 = []; b165 = []; c165 = [];

for j = 1:length(ngml)
	a121 = [a121 eval(strcat("Routs(j).",vals121,"(15+1)"))];
    a165 = [a165 eval(strcat("Routs(j+length(ngml)).",vals165,"(15+1)"))];
	b121 = [b121 eval(strcat("Routs(j).",vals121,"(60+1)"))];
    b165 = [b165 eval(strcat("Routs(j+length(ngml)).",vals165,"(60+1)"))];
	c121 = [c121 eval(strcat("Routs(j).",vals121,"(240+1)"))];
    c165 = [c165 eval(strcat("Routs(j+length(ngml)).",vals165,"(240+1)"))];
end

maxy = max(max([a121 a165 b121 b165 c121 c165]));
activeFig = figure('visible',visibility);
ax1 = subplot(1,3,1);
plot(ngml*v121ngmltonM,a121,'LineWidth',2.0); hold on;
plot(ngml*v165ngmltonM,a165,'LineWidth',2.0); 
legend('V121','V165','Location','northwest');
colororder(ligandmap);
title(ax1, '15min');
xlabel(ax1,'ligand (nM)');
set(ax1,'xscale','log');
xlim(ax1,[minngml maxngml]);
ylim(ax1,[0 maxy]);
ylabel(ax1,yaxistitle);
a = [ngml' a121' a165'];
csvname = strcat("outputDataDr2/",figname,"a.csv");
csvwrite(csvname,a);
  ax2 = subplot(1,3,2);
  plot(ngml*v121ngmltonM,b121,'LineWidth',2.0); hold on;
  plot(ngml*v165ngmltonM,b165,'LineWidth',2.0); 
  legend('V121','V165','Location','northwest');
  colororder(ligandmap);
  title(ax2, '60min');
  xlabel(ax2,'ligand (nM)');
  set(ax2,'xscale','log');
  xlim(ax2,[minngml maxngml]);
  ylim(ax2,[0 maxy]);
  ylabel(ax2,yaxistitle);
  b = [ngml' b121' b165'];
  csvname = strcat("outputDataDr2/",figname,"b.csv");
  csvwrite(csvname,b);
    ax3 = subplot(1,3,3);
    plot(ngml*v121ngmltonM,c121,'LineWidth',2.0); hold on;
    plot(ngml*v165ngmltonM,c165,'LineWidth',2.0); 
    legend('V121','V165','Location','northwest');
    colororder(ligandmap);
    title(ax3, '240min');
    xlabel(ax3,'ligand (nM)');
    set(ax3,'xscale','log');
    xlim(ax3,[minngml maxngml]);
    ylim(ax3,[0 maxy]);
    ylabel(ax3,yaxistitle);
    c = [ngml' c121' c165'];
    csvname = strcat("outputDataDr2/",figname,"c.csv");
    csvwrite(csvname,c);
set(activeFig, 'Position',[0 0 900 250])
set(activeFig, 'PaperUnits', 'inches');
set(activeFig, 'PaperSize', [4 6]);
exportgraphics(activeFig,strcat("outputFigsDr2/",figname,".png"),'Resolution',300);

end