function RDimznFig = VizRsdimerization(Routs,valsR1,valsR2,valsN1,figname,ngml,visibility,location,ligandmap)

minngml = min(ngml(2:end));
maxngml = max(ngml(2:end));

aR1 = []; bR1 = []; cR1 = []; dR1 = [];
aR2 = []; bR2 = []; cR2 = []; dR2 = [];
aN1 = []; bN1 = []; cN1 = []; dN1 = [];

for j = 1:length(ngml)
	aR1 = [aR1 eval(strcat("Routs(j).",valsR1,"(end)"))];
	aR2 = [aR2 eval(strcat("Routs(j).",valsR2,"(end)"))];
    aN1 = [aN1 eval(strcat("Routs(j).",valsN1,"(end)"))];
	bR1 = [bR1 eval(strcat("Routs(j+length(ngml)).",valsR1,"(end)"))];
	bR2 = [bR2 eval(strcat("Routs(j+length(ngml)).",valsR2,"(end)"))];
    bN1 = [bN1 eval(strcat("Routs(j+length(ngml)).",valsN1,"(end)"))];
	cR1 = [cR1 eval(strcat("Routs(j+2*length(ngml)).",valsR1,"(end)"))];
	cR2 = [cR2 eval(strcat("Routs(j+2*length(ngml)).",valsR2,"(end)"))];
    cN1 = [cN1 eval(strcat("Routs(j+2*length(ngml)).",valsN1,"(end)"))];
	dR1 = [dR1 eval(strcat("Routs(j+3*length(ngml)).",valsR1,"(end)"))];
	dR2 = [dR2 eval(strcat("Routs(j+3*length(ngml)).",valsR2,"(end)"))];
    dN1 = [dN1 eval(strcat("Routs(j+3*length(ngml)).",valsN1,"(end)"))];
end

RDimznFig = figure('visible',visibility);
ax1 = subplot(1,3,1);
plot(ngml,aR1','LineWidth',2.0); hold on;
plot(ngml,bR1','LineWidth',2.0); 
plot(ngml,cR1','LineWidth',2.0); 
plot(ngml,dR1','LineWidth',2.0); 
colororder(ligandmap);
legend('V121', 'V165','PlGF1','PlGF2','Location',location);
title(ax1, 'VEGFR1 - Dimerization');
xlabel(ax1,'ligand ng/ml');
ylabel(ax1,'Dimeric Fraction');
set(ax1,'xscale','log');
xlim(ax1,[minngml maxngml]);
ylim(ax1,[0 1]);
a = [ngml' aR1' bR1' cR1' dR1'];
csvname = strcat("outputDataDr2/",figname,"a.csv");
csvwrite(csvname,a);
  ax2 = subplot(1,3,2);
  plot(ngml,aR2','LineWidth',2.0); hold on;
  plot(ngml,bR2','LineWidth',2.0); 
  plot(ngml,cR2','LineWidth',2.0); 
  plot(ngml,dR2','LineWidth',2.0); 
  colororder(ligandmap);
  legend('V121', 'V165','PlGF1','PlGF2','Location',location);
  title(ax2, 'VEGFR2 - Dimerization');
  xlabel(ax2,'ligand ng/ml');
  ylabel(ax2,'Dimeric Fraction');
  set(ax2,'xscale','log');
  xlim(ax2,[minngml maxngml]);
  ylim(ax2,[0 1]);
  b = [ngml' aR2' bR2' cR2' dR2'];
  csvname = strcat("outputDataDr2/",figname,"b.csv");
  csvwrite(csvname,b);
    ax3 = subplot(1,3,3);
    plot(ngml,aN1','LineWidth',2.0); hold on;
    plot(ngml,bN1','LineWidth',2.0); 
    plot(ngml,cN1','LineWidth',2.0); 
    plot(ngml,dN1','LineWidth',2.0); 
    colororder(ligandmap);
    legend('V121', 'V165','PlGF1','PlGF2','Location',location);
    title(ax3, 'NRP1 - Dimerization');
    xlabel(ax3,'ligand ng/ml');
    ylabel(ax3,'Dimeric Fraction');
    set(ax3,'xscale','log');
    xlim(ax3,[minngml maxngml]);
    ylim(ax3,[0 1]);
    c = [ngml' aN1' bN1' cN1' dN1'];
    csvname = strcat("outputDataDr2/",figname,"c.csv");
    csvwrite(csvname,c);
set(RDimznFig, 'Position',[0 0 900 250])
set(RDimznFig, 'PaperUnits', 'inches');
set(RDimznFig, 'PaperSize', [4 6]);
exportgraphics(RDimznFig,strcat("outputFigsDr2/",figname,".png"),'Resolution',300);
end

