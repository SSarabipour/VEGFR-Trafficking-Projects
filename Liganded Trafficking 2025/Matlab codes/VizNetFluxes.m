
function NetFluxFig = VizNetFluxes(R1rates,receptorname,ngml,visibility)

minngml = min(ngml(2:end));
maxngml = max(ngml(2:end));

a0121 = []; a0165 = []; a0P2 = []; a0P1 = []; 
a15121 = []; a15165 = []; a15P2 = []; a15P1 = []; 
a60121 = []; a60165 = []; a60P2 = []; a60P1 = []; 
a240121 = []; a240165 = []; a240P2 = []; a240P1 = []; 
b0121 = []; b0165 = []; b0P2 = []; b0P1 = []; 
b15121 = []; b15165 = []; b15P2 = []; b15P1 = []; 
b60121 = []; b60165 = []; b60P2 = []; b60P1 = []; 
b240121 = []; b240165 = []; b240P2 = []; b240P1 = []; 
c0121 = []; c0165 = []; c0P2 = []; c0P1 = []; 
c15121 = []; c15165 = []; c15P2 = []; c15P1 = []; 
c60121 = []; c60165 = []; c60P2 = []; c60P1 = []; 
c240121 = []; c240165 = []; c240P2 = []; c240P1 = []; 


for j = 1:length(ngml)
	a0121 = [a0121 R1rates(j).Surf_net(0+1)];
    b0121 = [b0121 R1rates(j).Rab4_net(0+1)];
	c0121 = [c0121 R1rates(j).Rab11_net(0+1)];
	a0165 = [a0165 R1rates(j+length(ngml)).Surf_net(0+1)];
    b0165 = [b0165 R1rates(j+length(ngml)).Rab4_net(0+1)];
	c0165 = [c0165 R1rates(j+length(ngml)).Rab11_net(0+1)];
	a0P1 = [a0P1 R1rates(j+2*length(ngml)).Surf_net(0+1)];
    b0P1 = [b0P1 R1rates(j+2*length(ngml)).Rab4_net(0+1)];
	c0P1 = [c0P1 R1rates(j+2*length(ngml)).Rab11_net(0+1)];
	a0P2 = [a0P2 R1rates(j+3*length(ngml)).Surf_net(0+1)];
    b0P2 = [b0P2 R1rates(j+3*length(ngml)).Rab4_net(0+1)];
	c0P2 = [c0P2 R1rates(j+3*length(ngml)).Rab11_net(0+1)];

	a15121 = [a15121 R1rates(j).Surf_net(15+1)];
    b15121 = [b15121 R1rates(j).Rab4_net(15+1)];
	c15121 = [c15121 R1rates(j).Rab11_net(15+1)];
	a15165 = [a15165 R1rates(j+length(ngml)).Surf_net(15+1)];
    b15165 = [b15165 R1rates(j+length(ngml)).Rab4_net(15+1)];
	c15165 = [c15165 R1rates(j+length(ngml)).Rab11_net(15+1)];
	a15P1 = [a15P1 R1rates(j+2*length(ngml)).Surf_net(15+1)];
    b15P1 = [b15P1 R1rates(j+2*length(ngml)).Rab4_net(15+1)];
	c15P1 = [c15P1 R1rates(j+2*length(ngml)).Rab11_net(15+1)];
	a15P2 = [a15P2 R1rates(j+3*length(ngml)).Surf_net(15+1)];
    b15P2 = [b15P2 R1rates(j+3*length(ngml)).Rab4_net(15+1)];
	c15P2 = [c15P2 R1rates(j+3*length(ngml)).Rab11_net(15+1)];

	a60121 = [a60121 R1rates(j).Surf_net(60+1)];
    b60121 = [b60121 R1rates(j).Rab4_net(60+1)];
	c60121 = [c60121 R1rates(j).Rab11_net(60+1)];
	a60165 = [a60165 R1rates(j+length(ngml)).Surf_net(60+1)];
    b60165 = [b60165 R1rates(j+length(ngml)).Rab4_net(60+1)];
	c60165 = [c60165 R1rates(j+length(ngml)).Rab11_net(60+1)];
	a60P1 = [a60P1 R1rates(j+2*length(ngml)).Surf_net(60+1)];
    b60P1 = [b60P1 R1rates(j+2*length(ngml)).Rab4_net(60+1)];
	c60P1 = [c60P1 R1rates(j+2*length(ngml)).Rab11_net(60+1)];
	a60P2 = [a60P2 R1rates(j+3*length(ngml)).Surf_net(60+1)];
    b60P2 = [b60P2 R1rates(j+3*length(ngml)).Rab4_net(60+1)];
	c60P2 = [c60P2 R1rates(j+3*length(ngml)).Rab11_net(60+1)];

	a240121 = [a240121 R1rates(j).Surf_net(240+1)];
    b240121 = [b240121 R1rates(j).Rab4_net(240+1)];
	c240121 = [c240121 R1rates(j).Rab11_net(240+1)];
	a240165 = [a240165 R1rates(j+length(ngml)).Surf_net(240+1)];
    b240165 = [b240165 R1rates(j+length(ngml)).Rab4_net(240+1)];
	c240165 = [c240165 R1rates(j+length(ngml)).Rab11_net(240+1)];
	a240P1 = [a240P1 R1rates(j+2*length(ngml)).Surf_net(240+1)];
    b240P1 = [b240P1 R1rates(j+2*length(ngml)).Rab4_net(240+1)];
	c240P1 = [c240P1 R1rates(j+2*length(ngml)).Rab11_net(240+1)];
	a240P2 = [a240P2 R1rates(j+3*length(ngml)).Surf_net(240+1)];
    b240P2 = [b240P2 R1rates(j+3*length(ngml)).Rab4_net(240+1)];
	c240P2 = [c240P2 R1rates(j+3*length(ngml)).Rab11_net(240+1)];
    
end

miny = min(min([a0121 a0165 a0P2 a0P1 a15121 a15165 a15P2 a15P1 a60121 a60165 a60P2 a60P1 a240121 a240165 a240P2 a240P1; ...
b0121 b0165 b0P2 b0P1 b15121 b15165 b15P2 b15P1 b60121 b60165 b60P2 b60P1 b240121 b240165 b240P2 b240P1; ... 
c0121 c0165 c0P2 c0P1 c15121 c15165 c15P2 c15P1 c60121 c60165 c60P2 c60P1 c240121 c240165 c240P2 c240P1]));

maxy = max(max([a0121 a0165 a0P2 a0P1 a15121 a15165 a15P2 a15P1 a60121 a60165 a60P2 a60P1 a240121 a240165 a240P2 a240P1; ...
b0121 b0165 b0P2 b0P1 b15121 b15165 b15P2 b15P1 b60121 b60165 b60P2 b60P1 b240121 b240165 b240P2 b240P1; ... 
c0121 c0165 c0P2 c0P1 c15121 c15165 c15P2 c15P1 c60121 c60165 c60P2 c60P1 c240121 c240165 c240P2 c240P1]));

a = [ngml' a0121' b0121' c0121' a0165' b0165' c0165' a0P1' b0P1' c0P1' a0P2' b0P2' c0P2'];
csvname = strcat("outputDataDr2/NetRates_",receptorname,"_time0.csv");
csvwrite(csvname,a);
b = [ngml' a15121' b15121' c15121' a15165' b15165' c15165' a15P1' b15P1' c15P1' a15P2' b15P2' c15P2'];
csvname = strcat("outputDataDr2/NetRates_",receptorname,"_time15.csv");
csvwrite(csvname,b);
c = [ngml' a60121' b60121' c60121' a60165' b60165' c60165' a60P1' b60P1' c60P1' a60P2' b60P2' c60P2'];
csvname = strcat("outputDataDr2/NetRates_",receptorname,"_time60.csv");
csvwrite(csvname,c);
d = [ngml' a240121' b240121' c240121' a240165' b240165' c240165' a240P1' b240P1' c240P1' a240P2' b240P2' c240P2'];
csvname = strcat("outputDataDr2/NetRates_",receptorname,"_time240.csv");
csvwrite(csvname,d);

NetFluxFig = figure('visible',visibility);
ax1 = subplot(3,4,1);
plot(ngml,a15121,'LineWidth',2); 
hold on;
plot(ngml,b15121,'LineWidth',2); 
plot(ngml,c15121,'LineWidth',2); 
legend('surf','rab4','rab11','Location','northwest');
title(ax1, '15min - fluxes - V121');
xlabel(ax1,'ligand ng/ml');
set(ax1,'xscale','log');
ylabel(ax1,'Net Flux (in - out) (#/sec)');
xlim(ax1,[minngml maxngml]);
ylim(ax1,[miny maxy]);

ax2 = subplot(3,4,2);
plot(ngml,a15165,'LineWidth',2); 
hold on;
plot(ngml,b15165,'LineWidth',2); 
plot(ngml,c15165,'LineWidth',2); 
legend('surf','rab4','rab11','Location','northwest');
title(ax2, '15min - fluxes - V165');
xlabel(ax2,'ligand ng/ml');
set(ax2,'xscale','log');
ylabel(ax2,'Net Flux (in - out) (#/sec)');
xlim(ax2,[minngml maxngml]);
ylim(ax2,[miny maxy]);

ax3 = subplot(3,4,3);
plot(ngml,a15P1,'LineWidth',2); 
hold on;
plot(ngml,b15P1,'LineWidth',2); 
plot(ngml,c15P1,'LineWidth',2); 
legend('surf','rab4','rab11','Location','northwest');
title(ax3, '15min - fluxes - PlGF1');
xlabel(ax3,'ligand ng/ml');
set(ax3,'xscale','log');
ylabel(ax3,'Net Flux (in - out) (#/sec)');
xlim(ax3,[minngml maxngml]);
ylim(ax3,[miny maxy]);

ax4 = subplot(3,4,4);
plot(ngml,a15P2,'LineWidth',2); 
hold on;
plot(ngml,b15P2,'LineWidth',2); 
plot(ngml,c15P2,'LineWidth',2); 
legend('surf','rab4','rab11','Location','northwest');
title(ax4, '15min - fluxes - PlGF2');
xlabel(ax4,'ligand ng/ml');
set(ax4,'xscale','log');
ylabel(ax4,'Net Flux (in - out) (#/sec)');
xlim(ax4,[minngml maxngml]);
ylim(ax4,[miny maxy]);

ax5 = subplot(3,4,5);
plot(ngml,a60121,'LineWidth',2); 
hold on;
plot(ngml,b60121,'LineWidth',2); 
plot(ngml,c60121,'LineWidth',2); 
legend('surf','rab4','rab11','Location','northwest');
title(ax5, '60min - fluxes - V121');
xlabel(ax5,'ligand ng/ml');
set(ax5,'xscale','log');
ylabel(ax5,'Net Flux (in - out) (#/sec)');
xlim(ax5,[minngml maxngml]);
ylim(ax5,[miny maxy]);

ax6 = subplot(3,4,6);
plot(ngml,a60165,'LineWidth',2); 
hold on;
plot(ngml,b60165,'LineWidth',2); 
plot(ngml,c60165,'LineWidth',2); 
legend('surf','rab4','rab11','Location','northwest');
title(ax6, '60min - fluxes - V165');
xlabel(ax6,'ligand ng/ml');
set(ax6,'xscale','log');
ylabel(ax6,'Net Flux (in - out) (#/sec)');
xlim(ax6,[minngml maxngml]);
ylim(ax6,[miny maxy]);

ax7 = subplot(3,4,7) ;
plot(ngml,a60P1,'LineWidth',2); 
hold on;
plot(ngml,b60P1,'LineWidth',2); 
plot(ngml,c60P1,'LineWidth',2); 
legend('surf','rab4','rab11','Location','northwest');
title(ax7, '60min - fluxes - PlGF1');
xlabel(ax7,'ligand ng/ml');
set(ax7,'xscale','log');
ylabel(ax7,'Net Flux (in - out) (#/sec)');
xlim(ax7,[minngml maxngml]);
ylim(ax7,[miny maxy]);

ax8 = subplot(3,4,8);
plot(ngml,a60P2,'LineWidth',2); 
hold on;
plot(ngml,b60P2,'LineWidth',2); 
plot(ngml,c60P2,'LineWidth',2); 
legend('surf','rab4','rab11','Location','northwest');
title(ax8, '60min - fluxes - PlGF2');
xlabel(ax8,'ligand ng/ml');
set(ax8,'xscale','log');
ylabel(ax8,'Net Flux (in - out) (#/sec)');
xlim(ax8,[minngml maxngml]);
ylim(ax8,[miny maxy]);

ax9 = subplot(3,4,9);
plot(ngml,a240121,'LineWidth',2); 
hold on;
plot(ngml,b240121,'LineWidth',2); 
plot(ngml,c240121,'LineWidth',2); 
legend('surf','rab4','rab11','Location','northwest');
title(ax9, '240min - fluxes - V121');
xlabel(ax9,'ligand ng/ml');
set(ax9,'xscale','log');
ylabel(ax9,'Net Flux (in - out) (#/sec)');
xlim(ax9,[minngml maxngml]);
ylim(ax9,[miny maxy]);

ax10 = subplot(3,4,10);
plot(ngml,a240165,'LineWidth',2); 
hold on;
plot(ngml,b240165,'LineWidth',2); 
plot(ngml,c240165,'LineWidth',2); 
legend('surf','rab4','rab11','Location','northwest');
title(ax10, '240min - fluxes - V165');
xlabel(ax10,'ligand ng/ml');
set(ax10,'xscale','log');
ylabel(ax10,'Net Flux (in - out) (#/sec)');
xlim(ax10,[minngml maxngml]);
ylim(ax10,[miny maxy]);

ax11 = subplot(3,4,11);
plot(ngml,a240P1,'LineWidth',2); 
hold on;
plot(ngml,b240P1,'LineWidth',2); 
plot(ngml,c240P1,'LineWidth',2); 
legend('surf','rab4','rab11','Location','northwest');
title(ax11, '240min - fluxes - PlGF1');
xlabel(ax11,'ligand ng/ml');
set(ax11,'xscale','log');
ylabel(ax11,'Net Flux (in - out) (#/sec)');
xlim(ax11,[minngml maxngml]);
ylim(ax11,[miny maxy]);

ax12 = subplot(3,4,12);
plot(ngml,a240P2,'LineWidth',2); 
hold on;
plot(ngml,b240P2,'LineWidth',2); 
plot(ngml,c240P2,'LineWidth',2); 
legend('surf','rab4','rab11','Location','northwest');
title(ax12, '240min - fluxes - PlGF2');
xlabel(ax12,'ligand ng/ml');
set(ax12,'xscale','log');
ylabel(ax12,'Net Flux (in - out) (#/sec)');
xlim(ax12,[minngml maxngml]);
ylim(ax12,[miny maxy]);

set(NetFluxFig, 'Position',[0 0 1000 700])
set(NetFluxFig, 'PaperUnits', 'inches');
set(NetFluxFig, 'PaperSize', [4 6]);
exportgraphics(NetFluxFig,strcat("outputFigsDr2/NetRates_",receptorname,".png"),'Resolution',300);

end
