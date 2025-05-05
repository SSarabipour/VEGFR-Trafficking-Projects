function R1N1DimznFig = VizR1N1dimerization(Routs,valsR1,valsN1,figname,ngml,visibility,location,ligandmap)

minngml = min(ngml(2:end));
maxngml = max(ngml(2:end));

aR1 = []; bR1 = []; cR1 = []; dR1 = [];
aN1 = []; bN1 = []; cN1 = []; dN1 = [];

for j = 1:length(ngml)
	aR1 = [aR1 eval(strcat("Routs(j).",valsR1,"(end)"))];
    aN1 = [aN1 eval(strcat("Routs(j).",valsN1,"(end)"))];
	bR1 = [bR1 eval(strcat("Routs(j+length(ngml)).",valsR1,"(end)"))];
    bN1 = [bN1 eval(strcat("Routs(j+length(ngml)).",valsN1,"(end)"))];
	cR1 = [cR1 eval(strcat("Routs(j+2*length(ngml)).",valsR1,"(end)"))];
    cN1 = [cN1 eval(strcat("Routs(j+2*length(ngml)).",valsN1,"(end)"))];
	dR1 = [dR1 eval(strcat("Routs(j+3*length(ngml)).",valsR1,"(end)"))];
    dN1 = [dN1 eval(strcat("Routs(j+3*length(ngml)).",valsN1,"(end)"))];
end

R1N1DimznFig = figure('visible',visibility);
ax1 = subplot(1,2,1);
plot(ngml,aR1,'LineWidth',2.0); hold on;
plot(ngml,bR1,'LineWidth',2.0); 
plot(ngml,cR1,'LineWidth',2.0); 
plot(ngml,dR1,'LineWidth',2.0); 
colororder(ligandmap);
legend('V121', 'V165','PlGF1','PlGF2','Location',location);
title(ax1, 'R1N1 Dimerization - R1 in R1N1');
xlabel(ax1,'ligand ng/ml');
ylabel(ax1,'Fraction of R1 in R1N1');
set(ax1,'xscale','log');
xlim(ax1,[minngml maxngml]);
ylim(ax1,[0 1]);
a = [ngml' aR1' bR1' cR1' dR1'];
csvname = strcat("outputDataDr2/",figname,"a.csv");
csvwrite(csvname,a);
  ax2 = subplot(1,2,2);
  plot(ngml,aN1,'LineWidth',2.0); hold on;
  plot(ngml,bN1,'LineWidth',2.0); 
  plot(ngml,cN1,'LineWidth',2.0); 
  plot(ngml,dN1,'LineWidth',2.0); 
  legend('V121', 'V165','PlGF1','PlGF2','Location',location);
  colororder(ligandmap);
  title(ax2, 'R1N1 Dimerization - N1 in R1N1');
  xlabel(ax2,'ligand ng/ml');
  ylabel(ax2,'Fraction of N1 in R1N1');
  set(ax2,'xscale','log');
  xlim(ax2,[minngml maxngml]);
  ylim(ax2,[0 1]);
  b = [ngml' aN1' bN1' cN1' dN1'];
  csvname = strcat("outputDataDr2/",figname,"b.csv");
  csvwrite(csvname,b);
set(R1N1DimznFig, 'Position',[0 0 574 250])
set(R1N1DimznFig, 'PaperUnits', 'inches');
set(R1N1DimznFig, 'PaperSize', [4 6]);
exportgraphics(R1N1DimznFig,strcat("outputFigsDr2/",figname,".png"),'Resolution',300);

end

