function [R1_exp,R1_sim,R2_exp,R2_sim,N1_exp,N1_sim] = UnligandedHUVECExperimentsAll(params0,speciesInit,baselines,perturbs,recsurftargets,runModel,createFigs)

timestp_sim = 60*5; % 60 seconds
suppressligands = 1;

numexpts = 9;
R1_exp = zeros(numexpts,1);
R2_exp = zeros(numexpts,1);
N1_exp = zeros(numexpts,1);

R1_sim = zeros(numexpts,1);
R2_sim = zeros(numexpts,1);
N1_sim = zeros(numexpts,1);

% Initiate figures IFF figures are requested
if createFigs == 1
    CHX_by_output = figure;
    siRab4_by_output = figure;
    siRab11_by_output = figure;
    siRab4_11_by_output = figure;
    
	Surface_figure = figure;
    CHX_and_siRNA_overview = figure;
end

%% PART 1. QUANTIFIED WESTERNS following CHX treatment
% Expts #1,2,3,4,5 - 0.5,1,1.5,2,4 hr CHX
parameters=params0;

TimeLen = 20*3600 ;            % 20 hrs of no KD or ntRNA knockdown, no CHX
Timepoints  = 0:timestp_sim:TimeLen; 

[~, timepointsa, species_outa, observables_outa] = eval(append(runModel,"(Timepoints', speciesInit, parameters, 1)"));

speciesInit2 = species_outa(end,:); % starting point for CHX treatment

parameters.kR1prod = params0.kR1prod*perturbs.R1reductionByCHX;		%units are #/um^2.sec
parameters.kR2prod = params0.kR2prod*perturbs.R2reductionByCHX; 	%units are #/um^2.sec
parameters.kN1prod = params0.kN1prod*perturbs.N1reductionByCHX;  	%units are #/um^2.sec

TimeLen = 4*3600 ;               % 4 hour, no knockdown, + CHX - note, below 5 hrs, outputs are every 5 mins; above 5 hrs, it's every hour.
Timepoints  = 0:timestp_sim:TimeLen; 
[~, timepoints1, species_outb, observables_outb] = eval(append(runModel,"(Timepoints', speciesInit2, parameters, 1)"));

ConcOut1 = CalcOutputs(observables_outb,suppressligands);

i = 1;
j = 30*(60/timestp_sim)+1; % 30 min
R1_sim(i) = ConcOut1.R1_total(j);
R2_sim(i) = ConcOut1.R2_total(j);
N1_sim(i) = ConcOut1.N1_total(j);
R1_exp(i) = 0.80*baselines.R1_total;
R2_exp(i) = 0.87*baselines.R2_total;
N1_exp(i) = 0.95*baselines.N1_total;

i = 2;
j = 60*(60/timestp_sim)+1; % 60 min
R1_sim(i) = ConcOut1.R1_total(j);
R2_sim(i) = ConcOut1.R2_total(j);
N1_sim(i) = ConcOut1.N1_total(j);
R1_exp(i) = 0.35*baselines.R1_total;
R2_exp(i) = 0.60*baselines.R2_total;
N1_exp(i) = 1.11*baselines.N1_total;

i = 3;
j = 90*(60/timestp_sim)+1; % 90 min
R1_sim(i) = ConcOut1.R1_total(j);
R2_sim(i) = ConcOut1.R2_total(j);
N1_sim(i) = ConcOut1.N1_total(j);
R1_exp(i) = 0.29*baselines.R1_total;
R2_exp(i) = 0.59*baselines.R2_total;
N1_exp(i) = 1.04*baselines.N1_total;

i = 4;
j = 120*(60/timestp_sim)+1; % 120 min
R1_sim(i) = ConcOut1.R1_total(j);
R2_sim(i) = ConcOut1.R2_total(j);
N1_sim(i) = ConcOut1.N1_total(j);
R1_exp(i) = 0.20*baselines.R1_total;
R2_exp(i) = 0.40*baselines.R2_total;
N1_exp(i) = 1.05*baselines.N1_total;

i = 5;
j = 240*(60/timestp_sim)+1; % 240 min
R1_sim(i) = ConcOut1.R1_total(j);
R2_sim(i) = ConcOut1.R2_total(j);
N1_sim(i) = ConcOut1.N1_total(j);
R1_exp(i) = 0.15*baselines.R1_total;
R2_exp(i) = 0.09*baselines.R2_total;
N1_exp(i) = 1.01*baselines.N1_total;

timepoints1all = [timepointsa;timepoints1(2:end)+timepointsa(end)];
species_out1all = [species_outa;species_outb(2:end,:)];
observables_out1all = [observables_outa;observables_outb(2:end,:)];
ConcOut1all = CalcOutputs(observables_out1all,suppressligands);

exptimepoints1 = [30 60 90 120 240];


%% PART 2. QUANTIFIED WESTERNS following siRNA to Rab4a, Rab11a or both (knock down experiments)
% Expt #6-8 Rab4, Rab11, Rab4+11 siRNA treatment

% siRNA Rab4
parameters=params0;
parameters.kR1Rab4a = params0.kR1Rab4a*perturbs.Rab4reductionBysiRNA;
parameters.kR1N1Rab4a = params0.kR1N1Rab4a*perturbs.Rab4reductionBysiRNA;
parameters.kR2Rab4a = params0.kR2Rab4a*perturbs.Rab4reductionBysiRNA;
parameters.kN1Rab4a = params0.kN1Rab4a*perturbs.Rab4reductionBysiRNA;

TimeLen = 24*3600 ;            % 24 hrs of no KD or ntRNA knockdown
Timepoints  = 0:timestp_sim:TimeLen; 
[~, timepoints2, species_out, observables_out] = eval(append(runModel,"(Timepoints', speciesInit, parameters, 1)"));
% speciesInit2 = species_outa(end,:); % starting point for additional
%treatment if needed
ConcOut2 = CalcOutputs(observables_out,suppressligands);

i = 6;
R1_sim(i) = ConcOut2.R1_total(end);
R2_sim(i) = ConcOut2.R2_total(end);
N1_sim(i) = ConcOut2.N1_total(end);
R1_exp(i) = 1.0*baselines.R1_total;
R2_exp(i) = 1.0*baselines.R2_total;
N1_exp(i) = 1.0*baselines.N1_total;


% siRNA Rab11
parameters=params0;
parameters.kR1Rab4at11a = params0.kR1Rab4at11a*perturbs.Rab11reductionBysiRNA;
parameters.kR1N1Rab4at11a = params0.kR1N1Rab4at11a*perturbs.Rab11reductionBysiRNA;
parameters.kR2Rab4at11a = params0.kR2Rab4at11a*perturbs.Rab11reductionBysiRNA;
parameters.kN1Rab4at11a = params0.kN1Rab4at11a*perturbs.Rab11reductionBysiRNA;

TimeLen = 24*3600 ;            % 20 hrs of no KD or ntRNA knockdown, no CHX
Timepoints  = 0:timestp_sim:TimeLen; 
[~, timepoints3, species_out, observables_out] = eval(append(runModel,"(Timepoints', speciesInit, parameters, 1)"));
% speciesInit2 = species_outa(end,:); % starting point for additional
%treatment if needed
ConcOut3 = CalcOutputs(observables_out,suppressligands);

i = 7;
R1_sim(i) = ConcOut3.R1_total(end);
R2_sim(i) = ConcOut3.R2_total(end);
N1_sim(i) = ConcOut3.N1_total(end);
R1_exp(i) = 1.0*baselines.R1_total;
R2_exp(i) = 1.0*baselines.R2_total;
N1_exp(i) = 1.0*baselines.N1_total;


% siRNA Rab4+Rab11
parameters=params0;
parameters.kR1Rab4a = params0.kR1Rab4a*perturbs.Rab4reductionBysiRNA;
parameters.kR1N1Rab4a = params0.kR1N1Rab4a*perturbs.Rab4reductionBysiRNA;
parameters.kR2Rab4a = params0.kR2Rab4a*perturbs.Rab4reductionBysiRNA;
parameters.kN1Rab4a = params0.kN1Rab4a*perturbs.Rab4reductionBysiRNA;
parameters.kR1Rab4at11a = params0.kR1Rab4at11a*perturbs.Rab11reductionBysiRNA;
parameters.kR1N1Rab4at11a = params0.kR1N1Rab4at11a*perturbs.Rab11reductionBysiRNA;
parameters.kR2Rab4at11a = params0.kR2Rab4at11a*perturbs.Rab11reductionBysiRNA;
parameters.kN1Rab4at11a = params0.kN1Rab4at11a*perturbs.Rab11reductionBysiRNA;

TimeLen = 24*3600 ;            % 24 hrs of no KD or ntRNA knockdown, no CHX
Timepoints  = 0:timestp_sim:TimeLen; 
[~, timepoints4, species_out, observables_out] = eval(append(runModel,"(Timepoints', speciesInit, parameters, 1)"));
% speciesInit2 = species_outa(end,:); % starting point for additional
%treatment if needed
ConcOut4 = CalcOutputs(observables_out,suppressligands);

i = 8;
R1_sim(i) = ConcOut4.R1_total(end);
R2_sim(i) = ConcOut4.R2_total(end);
N1_sim(i) = ConcOut4.N1_total(end);
R1_exp(i) = 1.0*baselines.R1_total;
R2_exp(i) = 1.0*baselines.R2_total;
N1_exp(i) = 1.0*baselines.N1_total;

exptimepoints2 = [24*60]; exptimepoints3 = exptimepoints2; exptimepoints4 = exptimepoints2;


%% PART 3. INITIAL BALANCE OF SURF/INTERNAL-Biotinylation measurements
% Expt #9 no treatment

i = 9;

R1_exp(i) = 0.1;  % surface is about 10% of total R1
R2_exp(i) = 0.51; % surface is about 51% of total R2
N1_exp(i) = 0.74; % surface is about 74% of total N1
R1_sim(i) = baselines.R1_percSurf; % return the surface R1 only
R2_sim(i) = baselines.R2_percSurf; % return the surface R2 only
N1_sim(i) = baselines.N1_percSurf; % return the surface N1 only


%% Visualize results part 3
% Create figures IFF figures are requested
if createFigs == 1
    figure(Surface_figure);
    subplot(1,2,1);
    graphing = [recsurftargets(1) baselines.R1_surf recsurftargets(2) baselines.R2_surf recsurftargets(3) baselines.N1_surf];
    barcolors = [0 0 0; .5 .5 .5; 0 0 0; .5 .5 .5; 0 0 0; .5 .5 .5];
    b1 = bar(graphing);
    b1.FaceColor = 'flat';
    b1.CData = barcolors;
    set(gca,'XTickLabel',{'R1 Expt','R1 Sim','R2 Expt','R2 Sim','N1 Expt','N1 Sim'});
    title('Receptor levels on surface')
% add x-axis labels !

    subplot(1,2,2);
    graphing = [0.1 baselines.R1_percSurf 0.51 baselines.R2_percSurf 0.74 baselines.N1_percSurf];
    b2 = bar(graphing);
    b2.FaceColor = 'flat';
    b2.CData = barcolors;
    set(gca,'XTickLabel',{'R1 Expt','R1 Sim','R2 Expt','R2 Sim','N1 Expt','N1 Sim'});
    title('% of receptor on surface')

	saveas(Surface_figure,"Fig1_SurfRatios_Fit.png")

end


%% VISUALIZATION OF EXPERIMENTAL FITS
% Create figures IFF figures are requested
if createFigs == 1
figure(CHX_and_siRNA_overview);
    ax1=subplot(5,3,1);
        plot(timepoints1/3600,ConcOut1.R1_total);
        hold on;
        scatter(exptimepoints1/60,R1_exp(1:5));
        title('VEGFR1 wc - no siRNA', 'Interpreter','none');
    ax2=subplot(5,3,2);
        plot(timepoints1/3600,ConcOut1.R2_total);
        hold on;
        scatter(exptimepoints1/60,R2_exp(1:5));
        title('VEGFR2 wc - no siRNA', 'Interpreter','none');
    ax3=subplot(5,3,3);
        plot(timepoints1/3600,ConcOut1.N1_total);
        hold on;
        scatter(exptimepoints1/60,N1_exp(1:5));
        title('NRP1 wc - no siRNA', 'Interpreter','none');
    ax4=subplot(5,3,4);
        plot(timepoints2/3600,ConcOut2.R1_total);
        hold on;
       scatter(exptimepoints2/60,R1_exp(6));
        title('VEGFR1 wc - Rab4 siRNA', 'Interpreter','none');
    ax5=subplot(5,3,5);
        plot(timepoints2/3600,ConcOut2.R2_total);
        hold on;
       scatter(exptimepoints2/60,R2_exp(6));
        title('VEGFR2 wc - Rab4 siRNA', 'Interpreter','none');
    ax6=subplot(5,3,6);
        plot(timepoints2/3600,ConcOut2.N1_total);
        hold on;
        scatter(exptimepoints2/60,N1_exp(6));
        title('NRP1 wc - Rab4 siRNA', 'Interpreter','none');
    ax7=subplot(5,3,7);
        plot(timepoints3/3600,ConcOut3.R1_total);
        hold on;
        scatter(exptimepoints3/60,R1_exp(7));
        title('VEGFR1 wc - Rab11 siRNA', 'Interpreter','none');
    ax8=subplot(5,3,8);
        plot(timepoints3/3600,ConcOut3.R2_total);
        hold on;
        scatter(exptimepoints3/60,R2_exp(7));
        xlabel('time (hours) under CHX');
        title('VEGFR2 wc - Rab11 siRNA', 'Interpreter','none');
    ax9=subplot(5,3,9);
        plot(timepoints3/3600,ConcOut3.N1_total);
        hold on;
        scatter(exptimepoints3/60,N1_exp(7));
        title('NRP1 wc - Rab11 siRNA', 'Interpreter','none');
    ax10=subplot(5,3,10);
        plot(timepoints4/3600,ConcOut4.R1_total);
        hold on;
        scatter(exptimepoints4/60,R1_exp(8));
        title('VEGFR1 wc - Rab4+11 siRNA', 'Interpreter','none');
    ax11=subplot(5,3,11);
        plot(timepoints4/3600,ConcOut4.R2_total);
        hold on;
        scatter(exptimepoints4/60,R2_exp(8));
        xlabel('time (hours) under CHX');
        title('VEGFR2 wc - Rab4+11 siRNA', 'Interpreter','none');
    ax12=subplot(5,3,12);
        plot(timepoints4/3600,ConcOut4.N1_total);
        hold on;
        scatter(exptimepoints4/60,N1_exp(8));
        title('NRP1 wc - Rab4+11 siRNA', 'Interpreter','none');
    ax13=subplot(5,3,13);
        plot(timepoints1all/3600,ConcOut1all.R1_surf);
        hold on;
        plot(timepoints2/3600,ConcOut2.R1_surf);
        plot(timepoints3/3600,ConcOut3.R1_surf);
        plot(timepoints4/3600,ConcOut4.R1_surf);
        title('VEGFR1 surf - all conditions', 'Interpreter','none');
    ax14=subplot(5,3,14);
        plot(timepoints1all/3600,ConcOut1all.R2_surf);
        hold on;
        plot(timepoints2/3600,ConcOut2.R2_surf);
        plot(timepoints3/3600,ConcOut3.R2_surf);
        plot(timepoints4/3600,ConcOut4.R2_surf);
        title('VEGFR2 surf - all conditions', 'Interpreter','none');
    ax15=subplot(5,3,15);
        plot(timepoints1all/3600,ConcOut1all.N1_surf);
        hold on;
        plot(timepoints2/3600,ConcOut2.N1_surf);
        plot(timepoints3/3600,ConcOut3.N1_surf);
        plot(timepoints4/3600,ConcOut4.N1_surf);
        title('NRP1 surf - all conditions', 'Interpreter','none');% 	currentFigure = gcf;
%     title(currentFigure.Children(end),'Fit to Exptl Data');
    sgtitle('Fit to Exptl Data');
saveas(CHX_and_siRNA_overview,"Fig2_CHX_siRNA_Fit.png")
save('TimecourseCHXExpts.mat','timepoints1','ConcOut1','timepoints2','ConcOut2','timepoints3','ConcOut3',...
    'timepoints4','ConcOut4','exptimepoints1','R1_exp','R2_exp','N1_exp');

timecourses_forgraphing = [timepoints1/60 ConcOut1.R1_total'/max(ConcOut1.R1_total)*100 ...
                            ConcOut1.R2_total'/max(ConcOut1.R2_total)*100 ConcOut1.N1_total'/max(ConcOut1.N1_total)*100];
csvwrite('TimecourseCHXExpts.csv',timecourses_forgraphing);

locations_forgraphing = [baselines.R1_surf/baselines.R1_total; baselines.R1_Rab4/baselines.R1_total; baselines.R1_Rab11/baselines.R1_total; ...
                         baselines.R2_surf/baselines.R2_total; baselines.R2_Rab4/baselines.R2_total; baselines.R2_Rab11/baselines.R2_total; ...
                         baselines.N1_surf/baselines.N1_total; baselines.N1_Rab4/baselines.N1_total; baselines.N1_Rab11/baselines.N1_total];
csvwrite('Locations_baselines.csv',locations_forgraphing*100);

rabtotal_forgraphing = [baselines.R1_total/baselines.R1_total; ConcOut2.R1_total(end)/baselines.R1_total; ConcOut3.R1_total(end)/baselines.R1_total; ConcOut4.R1_total(end)/baselines.R1_total; ...
                        baselines.R2_total/baselines.R2_total; ConcOut2.R2_total(end)/baselines.R2_total; ConcOut3.R2_total(end)/baselines.R2_total; ConcOut4.R2_total(end)/baselines.R2_total; ...
                        baselines.N1_total/baselines.N1_total; ConcOut2.N1_total(end)/baselines.N1_total; ConcOut3.N1_total(end)/baselines.N1_total; ConcOut4.N1_total(end)/baselines.N1_total];
csvwrite('Totals_RabPerturbations.csv',rabtotal_forgraphing*100);

rabsurf_forgraphing = [baselines.R1_surf/baselines.R1_surf; ConcOut2.R1_surf(end)/baselines.R1_surf; ConcOut3.R1_surf(end)/baselines.R1_surf; ConcOut4.R1_surf(end)/baselines.R1_surf; ...
                       baselines.R2_surf/baselines.R2_surf; ConcOut2.R2_surf(end)/baselines.R2_surf; ConcOut3.R2_surf(end)/baselines.R2_surf; ConcOut4.R2_surf(end)/baselines.R2_surf; ...
                       baselines.N1_surf/baselines.N1_surf; ConcOut2.N1_surf(end)/baselines.N1_surf; ConcOut3.N1_surf(end)/baselines.N1_surf; ConcOut4.N1_surf(end)/baselines.N1_surf];
csvwrite('Surface_RabPerturbations.csv',rabsurf_forgraphing*100);

rabsurfperc_forgraphing = [baselines.R1_surf/baselines.R1_total; ConcOut2.R1_surf(end)/ConcOut2.R1_total(end); ConcOut3.R1_surf(end)/ConcOut3.R1_total(end); ConcOut4.R1_surf(end)/ConcOut4.R1_total(end); ...
                           baselines.R2_surf/baselines.R2_total; ConcOut2.R2_surf(end)/ConcOut2.R2_total(end); ConcOut3.R2_surf(end)/ConcOut3.R2_total(end); ConcOut4.R2_surf(end)/ConcOut4.R2_total(end); ...
                           baselines.N1_surf/baselines.N1_total; ConcOut2.N1_surf(end)/ConcOut2.N1_total(end); ConcOut3.N1_surf(end)/ConcOut3.N1_total(end); ConcOut4.N1_surf(end)/ConcOut4.N1_total(end)];
csvwrite('SurfacePerc_RabPerturbations.csv',rabsurfperc_forgraphing*100);

rabsurfperc_forgraphing_vs_ctrl = rabsurfperc_forgraphing;
rabsurfperc_forgraphing_vs_ctrl(1:4) = rabsurfperc_forgraphing_vs_ctrl(1:4)/rabsurfperc_forgraphing(1);
rabsurfperc_forgraphing_vs_ctrl(5:8) = rabsurfperc_forgraphing_vs_ctrl(5:8)/rabsurfperc_forgraphing(5);
rabsurfperc_forgraphing_vs_ctrl(9:12) = rabsurfperc_forgraphing_vs_ctrl(9:12)/rabsurfperc_forgraphing(9);
csvwrite('SurfacePercVsCtrl_RabPerturbations.csv',rabsurfperc_forgraphing_vs_ctrl*100);

% Visualize results parts 1-3

ConcOutNames = fieldnames(ConcOut1);
TotPanels = length(ConcOutNames);
PanelRows = ceil(length(ConcOutNames)/3);

figure(CHX_by_output);
for m=1:TotPanels
    eval(append("ax",num2str(m),"=subplot(",num2str(PanelRows),",3,",num2str(m),");"));
    plot(timepoints1/3600,ConcOut1.(ConcOutNames{m}));
    hold on;
    title(ConcOutNames{m}, 'Interpreter','none');
end
saveas(CHX_by_output,"Fig3A_CHXoutput.png")

figure(siRab4_by_output);
for m=1:TotPanels
    eval(append("ax",num2str(m),"=subplot(",num2str(PanelRows),",3,",num2str(m),");"));
    plot(timepoints2/3600,ConcOut2.(ConcOutNames{m}));
    hold on;
    title(ConcOutNames{m}, 'Interpreter','none');
end
saveas(siRab4_by_output,"Fig3B_siRab4_output.png")

figure(siRab11_by_output);
for m=1:TotPanels
    eval(append("ax",num2str(m),"=subplot(",num2str(PanelRows),",3,",num2str(m),");"));
    plot(timepoints3/3600,ConcOut3.(ConcOutNames{m}));
    hold on;
    title(ConcOutNames{m}, 'Interpreter','none');
end
saveas(siRab11_by_output,"Fig3C_siRab11_output.png")

figure(siRab4_11_by_output);
for m=1:TotPanels
    eval(append("ax",num2str(m),"=subplot(",num2str(PanelRows),",3,",num2str(m),");"));
    plot(timepoints4/3600,ConcOut4.(ConcOutNames{m}));
    hold on;
    title(ConcOutNames{m}, 'Interpreter','none');
end
saveas(siRab4_11_by_output,"Fig3D_siRab4_11_output.png")

end


end
