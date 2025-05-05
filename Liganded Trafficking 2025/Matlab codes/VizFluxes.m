
function FluxFig = VizFluxes(R1rates,receptorname,yaxistitle,ngml,visibility)

minngml = min(ngml(2:end));
maxngml = max(ngml(2:end));

a0121 = []; b0121 = []; c0121 = []; d0121 = []; e0121 = []; f0121 = [];
a0165 = []; b0165 = []; c0165 = []; d0165 = []; e0165 = []; f0165 = [];
a0P1 = []; b0P1 = []; c0P1 = []; d0P1 = []; e0P1 = []; f0P1 = [];
a0P2 = []; b0P2 = []; c0P2 = []; d0P2 = []; e0P2 = []; f0P2 = [];

a15121 = []; b15121 = []; c15121 = []; d15121 = []; e15121 = []; f15121 = [];
a15165 = []; b15165 = []; c15165 = []; d15165 = []; e15165 = []; f15165 = [];
a15P1 = []; b15P1 = []; c15P1 = []; d15P1 = []; e15P1 = []; f15P1 = [];
a15P2 = []; b15P2 = []; c15P2 = []; d15P2 = []; e15P2 = []; f15P2 = [];

a60121 = []; b60121 = []; c60121 = []; d60121 = []; e60121 = []; f60121 = [];
a60165 = []; b60165 = []; c60165 = []; d60165 = []; e60165 = []; f60165 = [];
a60P1 = []; b60P1 = []; c60P1 = []; d60P1 = []; e60P1 = []; f60P1 = [];
a60P2 = []; b60P2 = []; c60P2 = []; d60P2 = []; e60P2 = []; f60P2 = [];

a240121 = []; b240121 = []; c240121 = []; d240121 = []; e240121 = []; f240121 = [];
a240165 = []; b240165 = []; c240165 = []; d240165 = []; e240165 = []; f240165 = [];
a240P1 = []; b240P1 = []; c240P1 = []; d240P1 = []; e240P1 = []; f240P1 = [];
a240P2 = []; b240P2 = []; c240P2 = []; d240P2 = []; e240P2 = []; f240P2 = [];


for j = 1:length(ngml)
	a0121 = [a0121 R1rates(j).Surf_out_Rab4(0+1)];
    b0121 = [b0121 R1rates(j).Rab4_out_Rab7(0+1)];
	c0121 = [c0121 R1rates(j).Surf_in_Rab4(0+1)];
    d0121 = [d0121 R1rates(j).Rab11_in_Rab4(0+1)];
	e0121 = [e0121 R1rates(j).Surf_in_Rab11(0+1)];
    f0121 = [f0121 R1rates(j).Surf_in_Synth(0+1)];
    a0165 = [a0165 R1rates(j+length(ngml)).Surf_out_Rab4(0+1)];
    b0165 = [b0165 R1rates(j+length(ngml)).Rab4_out_Rab7(0+1)];
	c0165 = [c0165 R1rates(j+length(ngml)).Surf_in_Rab4(0+1)];
    d0165 = [d0165 R1rates(j+length(ngml)).Rab11_in_Rab4(0+1)];
	e0165 = [e0165 R1rates(j+length(ngml)).Surf_in_Rab11(0+1)];
    f0165 = [f0165 R1rates(j+length(ngml)).Surf_in_Synth(0+1)];
	a0P1 = [a0P1 R1rates(j+2*length(ngml)).Surf_out_Rab4(0+1)];
    b0P1 = [b0P1 R1rates(j+2*length(ngml)).Rab4_out_Rab7(0+1)];
	c0P1 = [c0P1 R1rates(j+2*length(ngml)).Surf_in_Rab4(0+1)];
    d0P1 = [d0P1 R1rates(j+2*length(ngml)).Rab11_in_Rab4(0+1)];
	e0P1 = [e0P1 R1rates(j+2*length(ngml)).Surf_in_Rab11(0+1)];
    f0P1 = [f0P1 R1rates(j+2*length(ngml)).Surf_in_Synth(0+1)];
	a0P2 = [a0P2 R1rates(j+3*length(ngml)).Surf_out_Rab4(0+1)];
    b0P2 = [b0P2 R1rates(j+3*length(ngml)).Rab4_out_Rab7(0+1)];
	c0P2 = [c0P2 R1rates(j+3*length(ngml)).Surf_in_Rab4(0+1)];
    d0P2 = [d0P2 R1rates(j+3*length(ngml)).Rab11_in_Rab4(0+1)];
	e0P2 = [e0P2 R1rates(j+3*length(ngml)).Surf_in_Rab11(0+1)];
    f0P2 = [f0P2 R1rates(j+3*length(ngml)).Surf_in_Synth(0+1)];

	a15121 = [a15121 R1rates(j).Surf_out_Rab4(15+1)];
    b15121 = [b15121 R1rates(j).Rab4_out_Rab7(15+1)];
	c15121 = [c15121 R1rates(j).Surf_in_Rab4(15+1)];
    d15121 = [d15121 R1rates(j).Rab11_in_Rab4(15+1)];
	e15121 = [e15121 R1rates(j).Surf_in_Rab11(15+1)];
    f15121 = [f15121 R1rates(j).Surf_in_Synth(15+1)];
    a15165 = [a15165 R1rates(j+length(ngml)).Surf_out_Rab4(15+1)];
    b15165 = [b15165 R1rates(j+length(ngml)).Rab4_out_Rab7(15+1)];
	c15165 = [c15165 R1rates(j+length(ngml)).Surf_in_Rab4(15+1)];
    d15165 = [d15165 R1rates(j+length(ngml)).Rab11_in_Rab4(15+1)];
	e15165 = [e15165 R1rates(j+length(ngml)).Surf_in_Rab11(15+1)];
    f15165 = [f15165 R1rates(j+length(ngml)).Surf_in_Synth(15+1)];
	a15P1 = [a15P1 R1rates(j+2*length(ngml)).Surf_out_Rab4(15+1)];
    b15P1 = [b15P1 R1rates(j+2*length(ngml)).Rab4_out_Rab7(15+1)];
	c15P1 = [c15P1 R1rates(j+2*length(ngml)).Surf_in_Rab4(15+1)];
    d15P1 = [d15P1 R1rates(j+2*length(ngml)).Rab11_in_Rab4(15+1)];
	e15P1 = [e15P1 R1rates(j+2*length(ngml)).Surf_in_Rab11(15+1)];
    f15P1 = [f15P1 R1rates(j+2*length(ngml)).Surf_in_Synth(15+1)];
	a15P2 = [a15P2 R1rates(j+3*length(ngml)).Surf_out_Rab4(15+1)];
    b15P2 = [b15P2 R1rates(j+3*length(ngml)).Rab4_out_Rab7(15+1)];
	c15P2 = [c15P2 R1rates(j+3*length(ngml)).Surf_in_Rab4(15+1)];
    d15P2 = [d15P2 R1rates(j+3*length(ngml)).Rab11_in_Rab4(15+1)];
	e15P2 = [e15P2 R1rates(j+3*length(ngml)).Surf_in_Rab11(15+1)];
    f15P2 = [f15P2 R1rates(j+3*length(ngml)).Surf_in_Synth(15+1)];

	a60121 = [a60121 R1rates(j).Surf_out_Rab4(60+1)];
    b60121 = [b60121 R1rates(j).Rab4_out_Rab7(60+1)];
	c60121 = [c60121 R1rates(j).Surf_in_Rab4(60+1)];
    d60121 = [d60121 R1rates(j).Rab11_in_Rab4(60+1)];
	e60121 = [e60121 R1rates(j).Surf_in_Rab11(60+1)];
    f60121 = [f60121 R1rates(j).Surf_in_Synth(60+1)];
    a60165 = [a60165 R1rates(j+length(ngml)).Surf_out_Rab4(60+1)];
    b60165 = [b60165 R1rates(j+length(ngml)).Rab4_out_Rab7(60+1)];
	c60165 = [c60165 R1rates(j+length(ngml)).Surf_in_Rab4(60+1)];
    d60165 = [d60165 R1rates(j+length(ngml)).Rab11_in_Rab4(60+1)];
	e60165 = [e60165 R1rates(j+length(ngml)).Surf_in_Rab11(60+1)];
    f60165 = [f60165 R1rates(j+length(ngml)).Surf_in_Synth(60+1)];
	a60P1 = [a60P1 R1rates(j+2*length(ngml)).Surf_out_Rab4(60+1)];
    b60P1 = [b60P1 R1rates(j+2*length(ngml)).Rab4_out_Rab7(60+1)];
	c60P1 = [c60P1 R1rates(j+2*length(ngml)).Surf_in_Rab4(60+1)];
    d60P1 = [d60P1 R1rates(j+2*length(ngml)).Rab11_in_Rab4(60+1)];
	e60P1 = [e60P1 R1rates(j+2*length(ngml)).Surf_in_Rab11(60+1)];
    f60P1 = [f60P1 R1rates(j+2*length(ngml)).Surf_in_Synth(60+1)];
	a60P2 = [a60P2 R1rates(j+3*length(ngml)).Surf_out_Rab4(60+1)];
    b60P2 = [b60P2 R1rates(j+3*length(ngml)).Rab4_out_Rab7(60+1)];
	c60P2 = [c60P2 R1rates(j+3*length(ngml)).Surf_in_Rab4(60+1)];
    d60P2 = [d60P2 R1rates(j+3*length(ngml)).Rab11_in_Rab4(60+1)];
	e60P2 = [e60P2 R1rates(j+3*length(ngml)).Surf_in_Rab11(60+1)];
    f60P2 = [f60P2 R1rates(j+3*length(ngml)).Surf_in_Synth(60+1)];

	a240121 = [a240121 R1rates(j).Surf_out_Rab4(240+1)];
    b240121 = [b240121 R1rates(j).Rab4_out_Rab7(240+1)];
	c240121 = [c240121 R1rates(j).Surf_in_Rab4(240+1)];
    d240121 = [d240121 R1rates(j).Rab11_in_Rab4(240+1)];
	e240121 = [e240121 R1rates(j).Surf_in_Rab11(240+1)];
    f240121 = [f240121 R1rates(j).Surf_in_Synth(240+1)];
    a240165 = [a240165 R1rates(j+length(ngml)).Surf_out_Rab4(240+1)];
    b240165 = [b240165 R1rates(j+length(ngml)).Rab4_out_Rab7(240+1)];
	c240165 = [c240165 R1rates(j+length(ngml)).Surf_in_Rab4(240+1)];
    d240165 = [d240165 R1rates(j+length(ngml)).Rab11_in_Rab4(240+1)];
	e240165 = [e240165 R1rates(j+length(ngml)).Surf_in_Rab11(240+1)];
    f240165 = [f240165 R1rates(j+length(ngml)).Surf_in_Synth(240+1)];
	a240P1 = [a240P1 R1rates(j+2*length(ngml)).Surf_out_Rab4(240+1)];
    b240P1 = [b240P1 R1rates(j+2*length(ngml)).Rab4_out_Rab7(240+1)];
	c240P1 = [c240P1 R1rates(j+2*length(ngml)).Surf_in_Rab4(240+1)];
    d240P1 = [d240P1 R1rates(j+2*length(ngml)).Rab11_in_Rab4(240+1)];
	e240P1 = [e240P1 R1rates(j+2*length(ngml)).Surf_in_Rab11(240+1)];
    f240P1 = [f240P1 R1rates(j+2*length(ngml)).Surf_in_Synth(240+1)];
	a240P2 = [a240P2 R1rates(j+3*length(ngml)).Surf_out_Rab4(240+1)];
    b240P2 = [b240P2 R1rates(j+3*length(ngml)).Rab4_out_Rab7(240+1)];
	c240P2 = [c240P2 R1rates(j+3*length(ngml)).Surf_in_Rab4(240+1)];
    d240P2 = [d240P2 R1rates(j+3*length(ngml)).Rab11_in_Rab4(240+1)];
	e240P2 = [e240P2 R1rates(j+3*length(ngml)).Surf_in_Rab11(240+1)];
    f240P2 = [f240P2 R1rates(j+3*length(ngml)).Surf_in_Synth(240+1)];

end

miny = min(min([a0121 b0121 c0121 d0121 e0121 f0121 a0165 b0165 c0165 d0165 e0165 f0165 a0P1 b0P1 c0P1 d0P1 e0P1 f0P1 a0P2 b0P2 c0P2 d0P2 e0P2 f0P2; ...
a15121 b15121 c15121 d15121 e15121 f15121 a15165 b15165 c15165 d15165 e15165 f15165 a15P1 b15P1 c15P1 d15P1 e15P1 f15P1 a15P2 b15P2 c15P2 d15P2 e15P2 f15P2; ...
a60121 b60121 c60121 d60121 e60121 f60121 a60165 b60165 c60165 d60165 e60165 f60165 a60P1 b60P1 c60P1 d60P1 e60P1 f60P1 a60P2 b60P2 c60P2 d60P2 e60P2 f60P2; ...
a240121 b240121 c240121 d240121 e240121 f240121 a240165 b240165 c240165 d240165 e240165 f240165 a240P1 b240P1 c240P1 d240P1 e240P1 f240P1 a240P2 b240P2 c240P2 d240P2 e240P2 f240P2]));

maxy = max(max([a0121 b0121 c0121 d0121 e0121 f0121 a0165 b0165 c0165 d0165 e0165 f0165 a0P1 b0P1 c0P1 d0P1 e0P1 f0P1 a0P2 b0P2 c0P2 d0P2 e0P2 f0P2; ...
a15121 b15121 c15121 d15121 e15121 f15121 a15165 b15165 c15165 d15165 e15165 f15165 a15P1 b15P1 c15P1 d15P1 e15P1 f15P1 a15P2 b15P2 c15P2 d15P2 e15P2 f15P2; ...
a60121 b60121 c60121 d60121 e60121 f60121 a60165 b60165 c60165 d60165 e60165 f60165 a60P1 b60P1 c60P1 d60P1 e60P1 f60P1 a60P2 b60P2 c60P2 d60P2 e60P2 f60P2; ...
a240121 b240121 c240121 d240121 e240121 f240121 a240165 b240165 c240165 d240165 e240165 f240165 a240P1 b240P1 c240P1 d240P1 e240P1 f240P1 a240P2 b240P2 c240P2 d240P2 e240P2 f240P2]));

a = [ngml' a0121' b0121' c0121' d0121' e0121' f0121' a0165' b0165' c0165' d0165' e0165' f0165' a0P1' b0P1' c0P1' d0P1' e0P1' f0P1' a0P2' b0P2' c0P2' d0P2' e0P2' f0P2'];
csvname = strcat("outputDataDr2/Rates_",receptorname,"_time0.csv");
csvwrite(csvname,a);
b = [ngml' a15121' b15121' c15121' d15121' e15121' f15121' a15165' b15165' c15165' d15165' e15165' f15165' a15P1' b15P1' c15P1' d15P1' e15P1' f15P1' a15P2' b15P2' c15P2' d15P2' e15P2' f15P2'];
csvname = strcat("outputDataDr2/Rates_",receptorname,"_time15.csv");
csvwrite(csvname,b);
c = [ngml' a60121' b60121' c60121' d60121' e60121' f60121' a60165' b60165' c60165' d60165' e60165' f60165' a60P1' b60P1' c60P1' d60P1' e60P1' f60P1' a60P2' b60P2' c60P2' d60P2' e60P2' f60P2'];
csvname = strcat("outputDataDr2/Rates_",receptorname,"_time60.csv");
csvwrite(csvname,c);
d = [ngml' a240121' b240121' c240121' d240121' e240121' f240121' a240165' b240165' c240165' d240165' e240165' f240165' a240P1' b240P1' c240P1' d240P1' e240P1' f240P1' a240P2' b240P2' c240P2' d240P2' e240P2' f240P2'];
csvname = strcat("outputDataDr2/Rates_",receptorname,"_time240.csv");
csvwrite(csvname,d);

FluxFig = figure('visible',visibility);
ax1 = subplot(3,4,1);
plot(ngml,a15121,'LineWidth',2); 
hold on;
plot(ngml,b15121,'LineWidth',2); 
plot(ngml,c15121,'LineWidth',2); 
plot(ngml,d15121,'LineWidth',2); 
plot(ngml,e15121,'LineWidth',2); 
plot(ngml,f15121,'LineWidth',2); 
legend('int','deg','rec4','tr411','rec11','prod','Location','northwest');
title(ax1, '15min - fluxes - V121');
xlabel(ax1,'ligand ng/ml');
set(ax1,'xscale','log');
ylabel(ax1,yaxistitle);
xlim(ax1,[minngml maxngml]);
ylim(ax1,[miny maxy]);

ax2 = subplot(3,4,2);
plot(ngml,a15165,'LineWidth',2); 
hold on;
plot(ngml,b15165,'LineWidth',2); 
plot(ngml,c15165,'LineWidth',2); 
plot(ngml,d15165,'LineWidth',2); 
plot(ngml,e15165,'LineWidth',2); 
plot(ngml,f15165,'LineWidth',2); 
legend('int','deg','rec4','tr411','rec11','prod','Location','northwest'); 
title(ax2, '15min - fluxes - V165');
xlabel(ax2,'ligand ng/ml');
set(ax2,'xscale','log');
ylabel(ax2,yaxistitle);
xlim(ax2,[minngml maxngml]);
ylim(ax2,[miny maxy]);

ax3 = subplot(3,4,3);
plot(ngml,a15P1,'LineWidth',2); 
hold on;
plot(ngml,b15P1,'LineWidth',2); 
plot(ngml,c15P1,'LineWidth',2); 
plot(ngml,d15P1,'LineWidth',2); 
plot(ngml,e15P1,'LineWidth',2); 
plot(ngml,f15P1,'LineWidth',2); 
legend('int','deg','rec4','tr411','rec11','prod','Location','northwest'); 
title(ax3, '15min - fluxes - PlGF1');
xlabel(ax3,'ligand ng/ml');
set(ax3,'xscale','log');
ylabel(ax3,yaxistitle);
xlim(ax3,[minngml maxngml]);
ylim(ax3,[miny maxy]);

ax4 = subplot(3,4,4);
plot(ngml,a15P2,'LineWidth',2); 
hold on;
plot(ngml,b15P2,'LineWidth',2); 
plot(ngml,c15P2,'LineWidth',2); 
plot(ngml,d15P2,'LineWidth',2); 
plot(ngml,e15P2,'LineWidth',2); 
plot(ngml,f15P2,'LineWidth',2); 
legend('int','deg','rec4','tr411','rec11','prod','Location','northwest'); 
title(ax4, '15min - fluxes - PlGF2');
xlabel(ax4,'ligand ng/ml');
set(ax4,'xscale','log');
ylabel(ax4,yaxistitle);
xlim(ax4,[minngml maxngml]);
ylim(ax4,[miny maxy]);

ax5 = subplot(3,4,5);
plot(ngml,a60121,'LineWidth',2); 
hold on;
plot(ngml,b60121,'LineWidth',2); 
plot(ngml,c60121,'LineWidth',2); 
plot(ngml,d60121,'LineWidth',2); 
plot(ngml,e60121,'LineWidth',2); 
plot(ngml,f60121,'LineWidth',2); 
legend('int','deg','rec4','tr411','rec11','prod','Location','northwest'); 
title(ax5, '60min - fluxes - V121');
xlabel(ax5,'ligand ng/ml');
set(ax5,'xscale','log');
ylabel(ax5,yaxistitle);
xlim(ax5,[minngml maxngml]);
ylim(ax5,[miny maxy]);

ax6 = subplot(3,4,6);
plot(ngml,a60165,'LineWidth',2); 
hold on;
plot(ngml,b60165,'LineWidth',2); 
plot(ngml,c60165,'LineWidth',2); 
plot(ngml,d60165,'LineWidth',2); 
plot(ngml,e60165,'LineWidth',2); 
plot(ngml,f60165,'LineWidth',2); 
legend('int','deg','rec4','tr411','rec11','prod','Location','northwest'); 
title(ax6, '60min - fluxes - V165');
xlabel(ax6,'ligand ng/ml');
set(ax6,'xscale','log');
ylabel(ax6,yaxistitle);
xlim(ax6,[minngml maxngml]);
ylim(ax6,[miny maxy]);

ax7 = subplot(3,4,7);
plot(ngml,a60P1,'LineWidth',2); 
hold on;
plot(ngml,b60P1,'LineWidth',2); 
plot(ngml,c60P1,'LineWidth',2); 
plot(ngml,d60P1,'LineWidth',2); 
plot(ngml,e60P1,'LineWidth',2); 
plot(ngml,f60P1,'LineWidth',2); 
legend('int','deg','rec4','tr411','rec11','prod','Location','northwest'); 
title(ax7, '60min - fluxes - PlGF1');
xlabel(ax7,'ligand ng/ml');
set(ax7,'xscale','log');
ylabel(ax7,yaxistitle);
xlim(ax7,[minngml maxngml]);
ylim(ax7,[miny maxy]);

ax8 = subplot(3,4,8);
plot(ngml,a60P2,'LineWidth',2); 
hold on;
plot(ngml,b60P2,'LineWidth',2); 
plot(ngml,c60P2,'LineWidth',2); 
plot(ngml,d60P2,'LineWidth',2); 
plot(ngml,e60P2,'LineWidth',2); 
plot(ngml,f60P2,'LineWidth',2); 
legend('int','deg','rec4','tr411','rec11','prod','Location','northwest'); 
title(ax8, '60min - fluxes - PlGF2');
xlabel(ax8,'ligand ng/ml');
set(ax8,'xscale','log');
ylabel(ax8,yaxistitle);
xlim(ax8,[minngml maxngml]);
ylim(ax8,[miny maxy]);

ax9 = subplot(3,4,9);
plot(ngml,a240121,'LineWidth',2); 
hold on;
plot(ngml,b240121,'LineWidth',2); 
plot(ngml,c240121,'LineWidth',2); 
plot(ngml,d240121,'LineWidth',2); 
plot(ngml,e240121,'LineWidth',2); 
plot(ngml,f240121,'LineWidth',2); 
legend('int','deg','rec4','tr411','rec11','prod','Location','northwest'); 
title(ax9, '240min - fluxes - V121');
xlabel(ax9,'ligand ng/ml');
set(ax9,'xscale','log');
ylabel(ax9,yaxistitle);
xlim(ax9,[minngml maxngml]);
ylim(ax9,[miny maxy]);

ax10 = subplot(3,4,10);
plot(ngml,a240165,'LineWidth',2); 
hold on;
plot(ngml,b240165,'LineWidth',2); 
plot(ngml,c240165,'LineWidth',2); 
plot(ngml,d240165,'LineWidth',2); 
plot(ngml,e240165,'LineWidth',2); 
plot(ngml,f240165,'LineWidth',2); 
legend('int','deg','rec4','tr411','rec11','prod','Location','northwest'); 
title(ax10, '240min - fluxes - V165');
xlabel(ax10,'ligand ng/ml');
set(ax10,'xscale','log');
ylabel(ax10,yaxistitle);
xlim(ax10,[minngml maxngml]);
ylim(ax10,[miny maxy]);

ax11 = subplot(3,4,11);
plot(ngml,a240P1,'LineWidth',2); 
hold on;
plot(ngml,b240P1,'LineWidth',2); 
plot(ngml,c240P1,'LineWidth',2); 
plot(ngml,d240P1,'LineWidth',2); 
plot(ngml,e240P1,'LineWidth',2); 
plot(ngml,f240P1,'LineWidth',2); 
legend('int','deg','rec4','tr411','rec11','prod','Location','northwest'); 
title(ax11, '240min - fluxes - PlGF1');
xlabel(ax11,'ligand ng/ml');
set(ax11,'xscale','log');
ylabel(ax11,yaxistitle);
xlim(ax11,[minngml maxngml]);
ylim(ax11,[miny maxy]);

ax12 = subplot(3,4,12);
plot(ngml,a240P2,'LineWidth',2); 
hold on;
plot(ngml,b240P2,'LineWidth',2); 
plot(ngml,c240P2,'LineWidth',2); 
plot(ngml,d240P2,'LineWidth',2); 
plot(ngml,e240P2,'LineWidth',2); 
plot(ngml,f240P2,'LineWidth',2); 
legend('int','deg','rec4','tr411','rec11','prod','Location','northwest'); 
title(ax12, '240min - fluxes - PlGF2');
xlabel(ax12,'ligand ng/ml');
set(ax12,'xscale','log');
ylabel(ax12,yaxistitle);
xlim(ax12,[minngml maxngml]);
ylim(ax12,[miny maxy]);

set(FluxFig, 'Position',[0 0 1000 700])
set(FluxFig, 'PaperUnits', 'inches');
set(FluxFig, 'PaperSize', [4 6]);
exportgraphics(FluxFig,strcat("outputFigsDr2/Rates_",receptorname,".png"),'Resolution',300);

end
