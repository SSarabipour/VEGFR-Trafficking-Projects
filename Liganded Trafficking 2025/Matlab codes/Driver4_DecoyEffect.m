close all;
clear all;

model = "All_Lig_Traff_model_Bionetgen_20231016"; %"TraffPhosph"
parameters = base_parameters(model);
% parameters.kR2prod  = 0; % TEST
% parameters.kR2Rab5a = 0; % TEST
% parameters.kR2N1Rab5a = 0; % TEST
parameters.kVegfR2Rab5a = parameters.kR2Rab5a * 3;
parameters.kVegfR2N1Rab5a = parameters.kVegfR2Rab5a;

ngml = [0 1 5 10 25 50 100 150];
v165a = [6843*1000*ngml./50]; % this is a vector of vegf165 concentrations. 
plgf1 = v165a/6843*10140; % this is a vector of plgf1 concentrations. Writing it this may makes sure the vectors are the same length
timepts = [0 15 30 60 120 240];

% vegflist = ["vegf 0ng/mL", "vegf 1ng/mL", "vegf 5ng/mL", "vegf 10ng/mL", "vegf 25ng/mL" ,"vegf 50ng/mL", "vegf 100ng/mL" ,"vegf 150ng/mL"];
% plgflist = ["plgf 0ng/mL", "plgf 1ng/mL", "plgf 5ng/mL", "plgf 10ng/mL", "plgf 25ng/mL" ,"plgf 50ng/mL", "plgf 100ng/mL" ,"plgf 150ng/mL"];
% timelist = ["0 min", "15 min", "30 min", "60 min", "120 min", "240 min"];

ReFitProductionRates = 1 ; % 1 for yes, 0 for no
R1surf_target = 1800;
R2surf_target = 4900; % TEST 9607;
N1surf_target = 68000;

visibility = 'off';

%% RUN NO-LIGAND STEADY STATES FOR THE FOUR SCENARIOS: R1/R2/N1; --/R2/N1; R1/--/N1; R1/R2/--

% STEADY STATE - ALL RECEPTORS PRESENT
disp('running st st R1/R2/N1')

if ReFitProductionRates == 1
    p2(1)=parameters.kR1prod;
    p2(2)=parameters.kR2prod;
    p2(3)=parameters.kN1prod;
    lb=p2*.001;     % lower bound
    ub=p2*1000;     % upper bound
    targets = [R1surf_target; R2surf_target; N1surf_target];
    options = optimoptions('lsqnonlin','Display','iter'); % display output
    optimalR = lsqnonlin(@CostFxn_ReceptorProduction,p2,lb,ub,options,targets,parameters,model);
    parameters.kR1prod = optimalR(1);
    parameters.kR2prod = optimalR(2);
    parameters.kN1prod = optimalR(3);
    optimalR
end

[timepoints, species_out, observables_out] = SimToSteadyState(parameters,model);
ConcOut=CalcOutputsL(observables_out);
ConcOutNames = fieldnames(ConcOut);
speciesInit = species_out(end,:);
fig_stst = VizConcOut(ConcOut,timepoints);

% STEADY STATE - no R1
disp('running st st --/R2/N1')
parameters_noR1 = parameters;
parameters_noR1.kR1prod = 0;
parameters_noR1.VEGFR1_0 = 0;

if ReFitProductionRates == 1
    p2(1)=parameters_noR1.kR1prod;
    p2(2)=parameters_noR1.kR2prod;
    p2(3)=parameters_noR1.kN1prod/1e3;
    lb=p2*.0001;     % lower bound
    ub=p2*10000;     % upper bound
    targets = [0; R2surf_target; N1surf_target];
    options = optimoptions('lsqnonlin','Display','iter'); % display output
    optimalR = lsqnonlin(@CostFxn_ReceptorProduction,p2,lb,ub,options,targets,parameters_noR1,model);
    parameters_noR1.kR2prod = optimalR(2);
    parameters_noR1.kN1prod = optimalR(3);
    optimalR
    optimalR(3)
end

[timepoints, species_out, observables_out] = SimToSteadyState(parameters_noR1,model);
ConcOut_noR1=CalcOutputsL(observables_out);
speciesInit_noR1 = species_out(end,:);
fig_stst_noR1 = VizConcOut(ConcOut_noR1,timepoints);

% STEADY STATE - no R2
disp('running st st R1/--/N1')
parameters_noR2 = parameters;
parameters_noR2.kR2prod = 0;
parameters_noR2.VEGFR2_0 = 0;

if ReFitProductionRates == 1
    p2(1)=parameters_noR2.kR1prod;
    p2(2)=parameters_noR2.kR2prod;
    p2(3)=parameters_noR2.kN1prod;
    lb=p2*.001;     % lower bound
    ub=p2*1000;     % upper bound
    targets = [R1surf_target; 0; N1surf_target];
    options = optimoptions('lsqnonlin','Display','iter'); % display output
    optimalR = lsqnonlin(@CostFxn_ReceptorProduction,p2,lb,ub,options,targets,parameters_noR2,model);
    parameters_noR2.kR1prod = optimalR(1);
    parameters_noR2.kN1prod = optimalR(3);
    optimalR
end

[timepoints, species_out, observables_out] = SimToSteadyState(parameters_noR2,model);
ConcOut_noR2=CalcOutputsL(observables_out);
speciesInit_noR2 = species_out(end,:);
fig_stst_noR2 = VizConcOut(ConcOut_noR2,timepoints);

% STEADY STATE - no N1
disp('running st st R1/R2/--')
parameters_noN1 = parameters;
parameters_noN1.kN1prod = 0;
parameters_noN1.NRP1_0 = 0;

if ReFitProductionRates == 1
    p2(1)=parameters_noN1.kR1prod;
    p2(2)=parameters_noN1.kR2prod;
    p2(3)=parameters_noN1.kN1prod;
    lb=p2*.001;     % lower bound
    ub=p2*1000;     % upper bound
    targets = [R1surf_target; R2surf_target; 0];
    options = optimoptions('lsqnonlin','Display','iter'); % display output
    optimalR = lsqnonlin(@CostFxn_ReceptorProduction,p2,lb,ub,options,targets,parameters_noN1,model);
    parameters_noN1.kR1prod = optimalR(1);
    parameters_noN1.kN2prod = optimalR(2);
    optimalR
end

[timepoints, species_out, observables_out] = SimToSteadyState(parameters_noN1,model);
ConcOut_noN1=CalcOutputsL(observables_out);
speciesInit_noN1 = species_out(end,:);
fig_stst_noN1 = VizConcOut(ConcOut_noN1,timepoints);


%% RUN "VEGF x PLGF" SCENARIOS FOR THE FOUR SCENARIOS: R1/R2/N1; --/R2/N1; R1/--/N1; R1/R2/--

disp('running VEGFxPlGF R1/R2/N1')
Routs = RunVEGFPlGFarray(model,speciesInit,parameters,'v165a','plgf1',v165a,plgf1,timepts);
disp('running VEGFxPlGF --/R2/N1')
Routs_noR1 = RunVEGFPlGFarray(model,speciesInit_noR1,parameters_noR1,'v165a','plgf1',v165a,plgf1,timepts);
disp('running VEGFxPlGF R1/--/N1')
Routs_noR2 = RunVEGFPlGFarray(model,speciesInit_noR2,parameters_noR2,'v165a','plgf1',v165a,plgf1,timepts);
disp('running VEGFxPlGF R1/R2/--')
Routs_noN1 = RunVEGFPlGFarray(model,speciesInit_noN1,parameters_noN1,'v165a','plgf1',v165a,plgf1,timepts);

% old format
% [R1_total, R2_total, N1_total, R1_surf, R2_surf, N1_surf, R1R1_165a_total, R2R2_165a_total, R1R1_165a_surf, R2R2_165a_surf] = RunVEGFPlGFarray(speciesInit,parameters,v165a,plgf,timepts);
% [R1_total_noR1, R2_total_noR1, N1_total_noR1, R1_surf_noR1, R2_surf_noR1, N1_surf_noR1, R1R1_165a_total_noR1, R2R2_165a_total_noR1, R1R1_165a_surf_noR1, R2R2_165a_surf_noR1] = RunVEGFPlGFarray(speciesInit_noR1,parameters_noR1,v165a,plgf,timepts);
% [R1_total_noR2, R2_total_noR2, N1_total_noR2, R1_surf_noR2, R2_surf_noR2, N1_surf_noR2, R1R1_165a_total_noR2, R2R2_165a_total_noR2, R1R1_165a_surf_noR2, R2R2_165a_surf_noR2] = RunVEGFPlGFarray(speciesInit_noR2,parameters_noR2,v165a,plgf,timepts);
% [R1_total_noN1, R2_total_noN1, N1_total_noN1, R1_surf_noN1, R2_surf_noN1, N1_surf_noN1, R1R1_165a_total_noN1, R2R2_165a_total_noN1, R1R1_165a_surf_noN1, R2R2_165a_surf_noN1] = RunVEGFPlGFarray(speciesInit_noN1,parameters_noN1,v165a,plgf,timepts);

outnames = fieldnames(Routs);

% two-color divergent map (down = purple; up = green)
mapr = [linspace(0.33,1,101)  linspace(.99,0,100)];
mapg = [linspace(0,1,101)   linspace(.99,.33,100)];
mapb = [linspace(0.33,1,101) linspace(0.99,0,100)];
map =[mapr' mapg' mapb'];

% one-color scale map (low = white; high = aqua blue)
mapr = linspace(1,66/255,201);
mapg = linspace(1,149/255,201); 
mapb = linspace(1,1,201); 
map1 =[mapr' mapg' mapb'];

%% OUTPUT HEATMAPS FOR THE METRICS - ABSOLUTE VALUES

disp('creating heatmaps - actual values')

% VEGF vs timepoints (for each deletion) (for each output)
for i = 1:length(outnames)
    a = squeeze(Routs.(outnames{i})(:,1,:));
    b = squeeze(Routs_noR1.(outnames{i})(:,1,:));
	c = squeeze(Routs_noR2.(outnames{i})(:,1,:));
	d = squeeze(Routs_noN1.(outnames{i})(:,1,:));
    a(isinf(a)) = NaN; b(isinf(b)) = NaN; c(isinf(c)) = NaN; d(isinf(d)) = NaN;
    scalemax = max(max([a b c d]));
    if scalemax == 0; scalemax = 1; end
    figname = strcat("vafig",num2str(i));
    eval(strcat(figname," = figure('visible',visibility);"));
    calcname = strcat("V_ROUT_fig_",num2str(i,'%02i'),"_",outnames{i});
    ax1 = subplot(4,1,1);
    h1 = heatmap(ngml,timepts,a'); % (:,1,:) vegf and timepoints; % (1,:,:) plgf and timepoints
    h1.Colormap = map1;
	h1.ColorLimits = [0,scalemax];
    h1.CellLabelFormat = '%.0f';
    h1.MissingDataColor = [0.9 0.9 0.9];
    h1.MissingDataLabel = "n/a";
    h1.FontSize = 14;
    h1.Title = 'Outputs - R1/R2/N1';
    h1.XLabel = 'VEGFconc (ng/ml)';
    h1.YLabel = 'Timepoint (min)';
	a = [0 ngml;timepts' a'];
    csvname = strcat("outputDataDr4/",calcname,"a.csv");
    csvwrite(csvname,a);
      ax2 = subplot(4,1,2);
      h2 = heatmap(ngml,timepts,b'); % (:,1,:) vegf and timepoints; % (1,:,:) plgf and timepoints
      h2.Colormap = map1;
  	  h2.ColorLimits = [0,scalemax];
      h2.CellLabelFormat = '%.0f';
      h2.MissingDataColor = [0.9 0.9 0.9];
      h2.MissingDataLabel = "n/a";
      h2.FontSize = 14;
      h2.Title = 'Outputs - --/R2/N1';
      h2.XLabel = 'VEGFconc (ng/ml)';
      h2.YLabel = 'Timepoint (min)';
      b = [0 ngml;timepts' b'];
      csvname = strcat("outputDataDr4/",calcname,"b.csv");
      csvwrite(csvname,b);
        ax3 = subplot(4,1,3);
        h3 = heatmap(ngml,timepts,c'); % (:,1,:) vegf and timepoints; % (1,:,:) plgf and timepoints
        h3.Colormap = map1;
        h3.ColorLimits = [0,scalemax];
        h3.CellLabelFormat = '%.0f';
        h3.MissingDataColor = [0.9 0.9 0.9];
        h3.MissingDataLabel = "n/a";
        h3.FontSize = 14;
        h3.Title = 'Outputs - R1/--/N1';
        h3.XLabel = 'VEGFconc (ng/ml)';
        h3.YLabel = 'Timepoint (min)';
        c = [0 ngml;timepts' c'];
        csvname = strcat("outputDataDr4/",calcname,"c.csv");
        csvwrite(csvname,c);
          ax4 = subplot(4,1,4);
          h4 = heatmap(ngml,timepts,d'); % (:,1,:) vegf and timepoints; % (1,:,:) plgf and timepoints
          h4.Colormap = map1;
          h4.ColorLimits = [0,scalemax];
          h4.CellLabelFormat = '%.0f';
          h4.MissingDataColor = [0.9 0.9 0.9];
          h4.MissingDataLabel = "n/a";
          h4.FontSize = 14;
          h4.Title = 'Outputs - R1/R2/--';
          h4.XLabel = 'VEGFconc (ng/ml)';
          h4.YLabel = 'Timepoint (min)';
          d = [0 ngml;timepts' d'];
          csvname = strcat("outputDataDr4/",calcname,"d.csv");
          csvwrite(csvname,d);
    set(eval(figname), 'Position',[0 0 1000 1000])
    set(eval(figname), 'PaperUnits', 'inches');
    set(eval(figname), 'PaperSize', [4 4]);
    exportgraphics(eval(figname),strcat("outputFigsDr4/",calcname,".png"),'Resolution',300);
end

% PlGF vs timepoints (for each deletion) (for each output)
for i = 1:length(outnames)
    a = squeeze(Routs.(outnames{i})(1,:,:));
    b = squeeze(Routs_noR1.(outnames{i})(1,:,:));
	c = squeeze(Routs_noR2.(outnames{i})(1,:,:));
	d = squeeze(Routs_noN1.(outnames{i})(1,:,:));
    a(isinf(a)) = NaN; b(isinf(b)) = NaN; c(isinf(c)) = NaN; d(isinf(d)) = NaN;
    scalemax = max(max([a b c d]));
    if scalemax == 0; scalemax = 1; end
    figname = strcat("pafig",num2str(i));
    eval(strcat(figname," = figure('visible',visibility);"));
    calcname = strcat("P_ROUT_fig_",num2str(i,'%02i'),"_",outnames{i});
    ax1 = subplot(4,1,1);
    h1 = heatmap(ngml,timepts,a'); % (:,1,:) vegf and timepoints; % (1,:,:) plgf and timepoints
    h1.Colormap = map1;
	h1.ColorLimits = [0,scalemax];
    h1.CellLabelFormat = '%.0f';
    h1.MissingDataColor = [0.9 0.9 0.9];
    h1.MissingDataLabel = "n/a";
    h1.FontSize = 14;
    h1.Title = 'Outputs - R1/R2/N1';
    h1.XLabel = 'PlGFconc (ng/ml)';
    h1.YLabel = 'Timepoint (min)';
	a = [0 ngml;timepts' a'];
    csvname = strcat("outputDataDr4/",calcname,"a.csv");
    csvwrite(csvname,a);
      ax2 = subplot(4,1,2);
      h2 = heatmap(ngml,timepts,b'); % (:,1,:) vegf and timepoints; % (1,:,:) plgf and timepoints
      h2.Colormap = map1;
	  h2.ColorLimits = [0,scalemax];
      h2.CellLabelFormat = '%.0f';
      h2.MissingDataColor = [0.9 0.9 0.9];
      h2.MissingDataLabel = "n/a";
      h2.FontSize = 14;
      h2.Title = 'Outputs - --/R2/N1';
      h2.XLabel = 'PlGFconc (ng/ml)';
      h2.YLabel = 'Timepoint (min)';
	  b = [0 ngml;timepts' b'];
      csvname = strcat("outputDataDr4/",calcname,"b.csv");
      csvwrite(csvname,b);
        ax3 = subplot(4,1,3);
        h3 = heatmap(ngml,timepts,c'); % (:,1,:) vegf and timepoints; % (1,:,:) plgf and timepoints
        h3.Colormap = map1;
        h3.ColorLimits = [0,scalemax];
        h3.CellLabelFormat = '%.0f';
        h3.MissingDataColor = [0.9 0.9 0.9];
        h3.MissingDataLabel = "n/a";
        h3.FontSize = 14;
        h3.Title = 'Outputs - R1/--/N1';
        h3.XLabel = 'PlGFconc (ng/ml)';
        h3.YLabel = 'Timepoint (min)';
        c = [0 ngml;timepts' c'];
        csvname = strcat("outputDataDr4/",calcname,"c.csv");
        csvwrite(csvname,c);
          ax4 = subplot(4,1,4);
          h4 = heatmap(ngml,timepts,d'); % (:,1,:) vegf and timepoints; % (1,:,:) plgf and timepoints
          h4.Colormap = map1;
          h4.ColorLimits = [0,scalemax];
          h4.CellLabelFormat = '%.0f';
          h4.MissingDataColor = [0.9 0.9 0.9];
          h4.MissingDataLabel = "n/a";
          h4.FontSize = 14;
          h4.Title = 'Outputs - R1/R2/--';
          h4.XLabel = 'PlGFconc (ng/ml)';
          h4.YLabel = 'Timepoint (min)';
          d = [0 ngml;timepts' d'];
          csvname = strcat("outputDataDr4/",calcname,"d.csv");
          csvwrite(csvname,d);
    set(eval(figname), 'Position',[0 0 1000 1000])
    set(eval(figname), 'PaperUnits', 'inches');
    set(eval(figname), 'PaperSize', [4 4]);
    exportgraphics(eval(figname),strcat("outputFigsDr4/",calcname,".png"),'Resolution',300);
end

%% CALCULATE THE DECOY EFFECTS - changes to outputs based on the

disp('calculating decoy effect metrics')
for i = 1:length(outnames)
    chR1.(outnames{i}) = Routs_noR1.(outnames{i}) - Routs.(outnames{i});
    relchR1.(outnames{i}) = chR1.(outnames{i}) ./ Routs.(outnames{i});
    relchR1.(outnames{i})(relchR1.(outnames{i})==Inf) = nan; 
    chR2.(outnames{i}) = Routs_noR2.(outnames{i}) - Routs.(outnames{i});
    relchR2.(outnames{i}) = chR2.(outnames{i}) ./ Routs.(outnames{i});
    relchR2.(outnames{i})(relchR2.(outnames{i})==Inf) = nan; 
    chN1.(outnames{i}) = Routs_noN1.(outnames{i}) - Routs.(outnames{i});
    relchN1.(outnames{i}) = chN1.(outnames{i}) ./ Routs.(outnames{i});
    relchN1.(outnames{i})(relchN1.(outnames{i})==Inf) = nan; 
end

%% OUTPUT HEATMAPS FOR THE DECOY EFFECT FOR VARIOUS OUTPUTS, VARIOUS TIMEPOINTS ETC

disp('creating heatmaps - decoy effect values')

% VEGF vs timepoints (for each deletion) (for each output)
for i = 1:length(outnames)
    a = squeeze(relchR1.(outnames{i})(:,1,:))*100;
	b = squeeze(relchR2.(outnames{i})(:,1,:))*100;
    c = squeeze(relchN1.(outnames{i})(:,1,:))*100;
    a(isinf(a)) = NaN; b(isinf(b)) = NaN; c(isinf(c)) = NaN; 
    figname = strcat("vfig",num2str(i));
    eval(strcat(figname," = figure('visible',visibility);"));
    calcname = strcat("Vt_fig_",num2str(i,'%02i'),"_",outnames{i});
    ax1 = subplot(3,1,1);
    h1 = heatmap(ngml,timepts,a','ColorLimits',[-100 100]); % (:,1,:) vegf and timepoints; % (1,:,:) plgf and timepoints
    h1.Colormap = map;
    h1.CellLabelFormat = '%.0f%%';
    h1.MissingDataColor = [0.9 0.9 0.9];
    h1.MissingDataLabel = "n/a";
    h1.FontSize = 14;
    h1.Title = 'Decoy Effect - effect of R1';
    h1.XLabel = 'VEGFconc (ng/ml)';
    h1.YLabel = 'Timepoint (min)';
    a = [0 ngml;timepts' a'];
    csvname = strcat("outputDataDr4/",calcname,"a.csv");
    csvwrite(csvname,a);
      ax2 = subplot(3,1,2);
      h2 = heatmap(ngml,timepts,b','ColorLimits',[-100 100]); % (:,1,:) vegf and timepoints; % (1,:,:) plgf and timepoints
      h2.Colormap = map;
      h2.CellLabelFormat = '%.0f%%';
      h2.MissingDataColor = [0.9 0.9 0.9];
      h2.MissingDataLabel = "n/a";
      h2.FontSize = 14;
      h2.Title = 'Decoy Effect - effect of R2';
      h2.XLabel = 'VEGFconc (ng/ml)';
      h2.YLabel = 'Timepoint (min)';
	  b = [0 ngml;timepts' b'];
      csvname = strcat("outputDataDr4/",calcname,"b.csv");
      csvwrite(csvname,b);
        ax3 = subplot(3,1,3);
        h3 = heatmap(ngml,timepts,c','ColorLimits',[-100 100]); % (:,1,:) vegf and timepoints; % (1,:,:) plgf and timepoints
        h3.Colormap = map;
        h3.CellLabelFormat = '%.0f%%';
        h3.MissingDataColor = [0.9 0.9 0.9];
        h3.MissingDataLabel = "n/a";
        h3.FontSize = 14;
        h3.Title = 'Decoy Effect - effect of N1';
        h3.XLabel = 'VEGFconc (ng/ml)';
        h3.YLabel = 'Timepoint (min)';
        c = [0 ngml;timepts' c'];
        csvname = strcat("outputDataDr4/",calcname,"c.csv");
        csvwrite(csvname,c);
    set(eval(figname), 'Position',[0 0 1000 729])
    set(eval(figname), 'PaperUnits', 'inches');
    set(eval(figname), 'PaperSize', [4 4]);
    exportgraphics(eval(figname),strcat("outputFigsDr4/",calcname,".png"),'Resolution',300);
end

% PlGF vs timepoints (for each deletion) (for each output)
for i = 1:length(outnames)
    a = squeeze(relchR1.(outnames{i})(1,:,:))*100;
    b = squeeze(relchR2.(outnames{i})(1,:,:))*100;
	c = squeeze(relchN1.(outnames{i})(1,:,:))*100;
    a(isinf(a)) = NaN; b(isinf(b)) = NaN; c(isinf(c)) = NaN; 
    figname = strcat("pfig",num2str(i));
    eval(strcat(figname," = figure('visible',visibility);"));
    calcname = strcat("Pt_fig_",num2str(i,'%02i'),"_",outnames{i});
    ax1 = subplot(3,1,1);
    h1 = heatmap(ngml,timepts,a','ColorLimits',[-100 100]); % (:,1,:) vegf and timepoints; % (1,:,:) plgf and timepoints
    h1.Colormap = map;
    h1.CellLabelFormat = '%.0f%%';
    h1.MissingDataColor = [0.9 0.9 0.9];
    h1.MissingDataLabel = "n/a";
    h1.FontSize = 14;
    h1.Title = 'Decoy Effect - effect of R1';
    h1.XLabel = 'PlGFconc (ng/ml)';
    h1.YLabel = 'Timepoint (min)';
    a = [0 ngml;timepts' a'];
    csvname = strcat("outputDataDr4/",calcname,"a.csv");
    csvwrite(csvname,a);
      ax2 = subplot(3,1,2);
      h2 = heatmap(ngml,timepts,b','ColorLimits',[-100 100]); % (:,1,:) vegf and timepoints; % (1,:,:) plgf and timepoints
      h2.Colormap = map;
      h2.CellLabelFormat = '%.0f%%';
      h2.MissingDataColor = [0.9 0.9 0.9];
      h2.MissingDataLabel = "n/a";
      h2.FontSize = 14;
      h2.Title = 'Decoy Effect - effect of R2';
      h2.XLabel = 'PlGFconc (ng/ml)';
      h2.YLabel = 'Timepoint (min)';
      b = [0 ngml;timepts' b'];
      csvname = strcat("outputDataDr4/",calcname,"b.csv");
      csvwrite(csvname,b);
        ax3 = subplot(3,1,3);
        h3 = heatmap(ngml,timepts,c','ColorLimits',[-100 100]); % (:,1,:) vegf and timepoints; % (1,:,:) plgf and timepoints
        h3.Colormap = map;
        h3.CellLabelFormat = '%.0f%%';
        h3.MissingDataColor = [0.9 0.9 0.9];
        h3.MissingDataLabel = "n/a";
        h3.FontSize = 14;
        h3.Title = 'Decoy Effect - effect of N1';
        h3.XLabel = 'PlGFconc (ng/ml)';
        h3.YLabel = 'Timepoint (min)';
        c = [0 ngml;timepts' c'];
        csvname = strcat("outputDataDr4/",calcname,"c.csv");
        csvwrite(csvname,c);
    set(eval(figname), 'Position',[0 0 1000 729])
    set(eval(figname), 'PaperUnits', 'inches');
    set(eval(figname), 'PaperSize', [4 4]);
    exportgraphics(eval(figname),strcat("outputFigsDr4/",calcname,".png"),'Resolution',300);
end

% VEGF vs outputs (for each timepoint)
for i = 1:length(timepts)
    a = zeros(length(v165a),length(outnames));
    for j = 1:length(outnames)
        a(:,j) = squeeze(relchR1.(outnames{j})(:,1,i))*100;
    end
    b = zeros(length(v165a),length(outnames));
    for j = 1:length(outnames)
        b(:,j) = squeeze(relchR2.(outnames{j})(:,1,i))*100;
    end
    c = zeros(length(v165a),length(outnames));
    for j = 1:length(outnames)
        c(:,j) = squeeze(relchN1.(outnames{j})(:,1,i))*100;
    end
    a(isinf(a)) = NaN; b(isinf(b)) = NaN; c(isinf(c)) = NaN; 
    figname = strcat("vofig",num2str(i));
    eval(strcat(figname," = figure('visible',visibility);"));
    calcname = strcat("Vo_fig_",num2str(i,'%02i'),"_",num2str(timepts(i),'%03i'),"mins");
    ax1 = subplot(1,3,1);
    h1 = heatmap(ngml,outnames,a','ColorLimits',[-100 100]); % (:,1,:) vegf and timepoints; % (1,:,:) plgf and timepoints
    h1.Colormap = map;
    h1.CellLabelFormat = '%.0f%%';
    h1.MissingDataColor = [0.9 0.9 0.9];
    h1.MissingDataLabel = "n/a";
    h1.FontSize = 14;
    h1.Title = 'Decoy Effect - effect of R1';
    h1.XLabel = 'VEGFconc (ng/ml)';
    h1.YLabel = 'Outputs';
    a = [ngml;a'];
    csvname = strcat("outputDataDr4/",calcname,"a.csv");
    csvwrite(csvname,a);
      ax2 = subplot(1,3,2);
      h2 = heatmap(ngml,outnames,b','ColorLimits',[-100 100]); % (:,1,:) vegf and timepoints; % (1,:,:) plgf and timepoints
      h2.Colormap = map;
      h2.CellLabelFormat = '%.0f%%';
      h2.MissingDataColor = [0.9 0.9 0.9];
      h2.MissingDataLabel = "n/a";
      h2.FontSize = 14;
      h2.Title = 'Decoy Effect - effect of R2';
      h2.XLabel = 'VEGFconc (ng/ml)';
      h2.YLabel = 'Outputs';
      b = [ngml;b'];
      csvname = strcat("outputDataDr4/",calcname,"b.csv");
      csvwrite(csvname,b);
        ax3 = subplot(1,3,3);
        h3 = heatmap(ngml,outnames,c','ColorLimits',[-100 100]); % (:,1,:) vegf and timepoints; % (1,:,:) plgf and timepoints
        h3.Colormap = map;
        h3.CellLabelFormat = '%.0f%%';
        h3.MissingDataColor = [0.9 0.9 0.9];
        h3.MissingDataLabel = "n/a";
        h3.FontSize = 14;
        h3.Title = 'Decoy Effect - effect of N1';
        h3.XLabel = 'VEGFconc (ng/ml)';
        h3.YLabel = 'Outputs';
        c = [ngml;c'];
        csvname = strcat("outputDataDr4/",calcname,"c.csv");
        csvwrite(csvname,c);
    set(eval(figname), 'Position',[0 0 1500 1000])
    set(eval(figname), 'PaperUnits', 'inches');
    set(eval(figname), 'PaperSize', [4 4]);
    exportgraphics(eval(figname),strcat("outputFigsDr4/",calcname,".png"),'Resolution',300);
end

% PlGF vs outputs (for each timepoint)
for i = 1:length(timepts)
    a = zeros(length(plgf1),length(outnames));
    for j = 1:length(outnames)
        a(:,j) = squeeze(relchR1.(outnames{j})(1,:,i))*100;
    end
    b = zeros(length(plgf1),length(outnames));
    for j = 1:length(outnames)
        b(:,j) = squeeze(relchR2.(outnames{j})(1,:,i))*100;
    end
    c = zeros(length(plgf1),length(outnames));
    for j = 1:length(outnames)
        c(:,j) = squeeze(relchN1.(outnames{j})(1,:,i))*100;
    end
    a(isinf(a)) = NaN; b(isinf(b)) = NaN; c(isinf(c)) = NaN; 
    figname = strcat("pofig",num2str(i));
    eval(strcat(figname," = figure('visible',visibility);"));
    calcname = strcat("Po_fig_",num2str(i,'%02i'),"_",num2str(timepts(i),'%03i'),"mins");
    ax1 = subplot(1,3,1);
    h1 = heatmap(ngml,outnames,a','ColorLimits',[-100 100]); % (:,1,:) vegf and timepoints; % (1,:,:) plgf and timepoints
    h1.Colormap = map;
    h1.CellLabelFormat = '%.0f%%';
    h1.MissingDataColor = [0.9 0.9 0.9];
    h1.MissingDataLabel = "n/a";
    h1.FontSize = 14;
    h1.Title = 'Decoy Effect - effect of R1';
    h1.XLabel = 'PlGFconc (ng/ml)';
    h1.YLabel = 'Outputs';
    a = [ngml;a'];
    csvname = strcat("outputDataDr4/",calcname,"a.csv");
    csvwrite(csvname,a);
      ax2 = subplot(1,3,2);
      h2 = heatmap(ngml,outnames,b','ColorLimits',[-100 100]); % (:,1,:) vegf and timepoints; % (1,:,:) plgf and timepoints
      h2.Colormap = map;
      h2.CellLabelFormat = '%.0f%%';
      h2.MissingDataColor = [0.9 0.9 0.9];
      h2.MissingDataLabel = "n/a";
      h2.FontSize = 14;
      h2.Title = 'Decoy Effect - effect of R2';
      h2.XLabel = 'PlGFconc (ng/ml)';
      h2.YLabel = 'Outputs';
      b = [ngml;b'];
      csvname = strcat("outputDataDr4/",calcname,"b.csv");
      csvwrite(csvname,b);
        ax3 = subplot(1,3,3);
        h3 = heatmap(ngml,outnames,c','ColorLimits',[-100 100]); % (:,1,:) vegf and timepoints; % (1,:,:) plgf and timepoints
        h3.Colormap = map;
        h3.CellLabelFormat = '%.0f%%';
        h3.MissingDataColor = [0.9 0.9 0.9];
        h3.MissingDataLabel = "n/a";
        h3.FontSize = 14;
        h3.Title = 'Decoy Effect - effect of N1';
        h3.XLabel = 'PlGFconc (ng/ml)';
        h3.YLabel = 'Outputs';
        c = [ngml;c'];
        csvname = strcat("outputDataDr4/",calcname,"c.csv");
        csvwrite(csvname,c);
    set(eval(figname), 'Position',[0 0 1500 1000])
    set(eval(figname), 'PaperUnits', 'inches');
    set(eval(figname), 'PaperSize', [4 4]);
    exportgraphics(eval(figname),strcat("outputFigsDr4/",calcname,".png"),'Resolution',300);
end

% VEGF vs PlGF (for each deletion) (for each output) - first 15 mins
tsel = 2; %[timepts(2) = 15]

for i = 1:length(outnames)
    a = squeeze(Routs.(outnames{i})(:,:,tsel));
    b = squeeze(Routs_noR1.(outnames{i})(:,:,tsel));
    c = squeeze(Routs_noR2.(outnames{i})(:,:,tsel));
	d = squeeze(Routs_noN1.(outnames{i})(:,:,tsel));
    a(isinf(a)) = NaN; b(isinf(b)) = NaN; c(isinf(c)) = NaN; d(isinf(d)) = NaN;
    scalemax = max(max([a b c d]));
    if scalemax == 0; scalemax = 1; end
    figname = strcat("vpfig",num2str(i));
    eval(strcat(figname," = figure('visible',visibility);"));
    calcname = strcat("VP_15m_ROUT_fig_",num2str(i,'%02i'),"_",outnames{i});
    ax1 = subplot(4,1,1);
    h1 = heatmap(ngml,ngml,a'); % (:,1,:) vegf and timepoints; % (1,:,:) plgf and timepoints
    h1.Colormap = map1;
	h1.ColorLimits = [0,scalemax];
	h1.CellLabelFormat = '%.0f';
    h1.MissingDataColor = [0.9 0.9 0.9];
    h1.MissingDataLabel = "n/a";
    h1.FontSize = 14;
    h1.Title = 'Outputs - R1/R2/N1';
    h1.XLabel = 'VEGFconc (ng/ml)';
    h1.YLabel = 'PlGFconc (ng/ml)';
    a = [0 ngml;ngml' a'];
    csvname = strcat("outputDataDr4/",calcname,"a.csv");
    csvwrite(csvname,a);
      ax2 = subplot(4,1,2);
      h2 = heatmap(ngml,ngml,b'); % (:,1,:) vegf and timepoints; % (1,:,:) plgf and timepoints
      h2.Colormap = map1;
      h2.ColorLimits = [0,scalemax];
      h2.CellLabelFormat = '%.0f';
      h2.MissingDataColor = [0.9 0.9 0.9];
      h2.MissingDataLabel = "n/a";
      h2.FontSize = 14;
      h2.Title = 'Outputs - --/R2/N1';
      h2.XLabel = 'VEGFconc (ng/ml)';
      h2.YLabel = 'PlGFconc (ng/ml)';
      b = [0 ngml;ngml' b'];
      csvname = strcat("outputDataDr4/",calcname,"b.csv");
      csvwrite(csvname,b);
        ax3 = subplot(4,1,3);
        h3 = heatmap(ngml,ngml,c'); % (:,1,:) vegf and timepoints; % (1,:,:) plgf and timepoints
        h3.Colormap = map1;
        h3.ColorLimits = [0,scalemax];
    	h3.CellLabelFormat = '%.0f';
        h3.MissingDataColor = [0.9 0.9 0.9];
        h3.MissingDataLabel = "n/a";
        h3.FontSize = 14;
        h3.Title = 'Outputs - R1/--/N1';
        h3.XLabel = 'VEGFconc (ng/ml)';
        h3.YLabel = 'PlGFconc (ng/ml)';
        c = [0 ngml;ngml' c'];
        csvname = strcat("outputDataDr4/",calcname,"c.csv");
        csvwrite(csvname,c);
          ax4 = subplot(4,1,4);
          h4 = heatmap(ngml,ngml,d'); % (:,1,:) vegf and timepoints; % (1,:,:) plgf and timepoints
          h4.Colormap = map1;
          h4.ColorLimits = [0,scalemax];
      	  h4.CellLabelFormat = '%.0f';
          h4.MissingDataColor = [0.9 0.9 0.9];
          h4.MissingDataLabel = "n/a";
          h4.FontSize = 14;
          h4.Title = 'Outputs - R1/R2/--';
          h4.XLabel = 'VEGFconc (ng/ml)';
          h4.YLabel = 'PlGFconc (ng/ml)';
          d = [0 ngml;ngml' d'];
          csvname = strcat("outputDataDr4/",calcname,"d.csv");
          csvwrite(csvname,d);
    set(eval(figname), 'Position',[0 0 1000 1000])
    set(eval(figname), 'PaperUnits', 'inches');
    set(eval(figname), 'PaperSize', [4 4]);
    exportgraphics(eval(figname),strcat("outputFigsDr4/",calcname,".png"),'Resolution',300);
end

% VEGF vs PlGF (for each deletion) (for each output) - first 15 mins
tsel = 5; %[timepts(5) = 120]

for i = 1:length(outnames)
    a = squeeze(Routs.(outnames{i})(:,:,tsel));
	b = squeeze(Routs_noR1.(outnames{i})(:,:,tsel));
	c = squeeze(Routs_noR2.(outnames{i})(:,:,tsel));
	d = squeeze(Routs_noN1.(outnames{i})(:,:,tsel));
    a(isinf(a)) = NaN; b(isinf(b)) = NaN; c(isinf(c)) = NaN; d(isinf(d)) = NaN;
    scalemax = max(max([a b c d]));
    if scalemax == 0; scalemax = 1; end
    figname = strcat("vp2fig",num2str(i));
    eval(strcat(figname," = figure('visible',visibility);"));
    calcname = strcat("VP_120m_ROUT_fig_",num2str(i,'%02i'),"_",outnames{i});
    ax1 = subplot(4,1,1);
    h1 = heatmap(ngml,ngml,a'); % (:,1,:) vegf and timepoints; % (1,:,:) plgf and timepoints
    h1.Colormap = map1;
	h1.ColorLimits = [0,scalemax];
	h1.CellLabelFormat = '%.0f';
    h1.MissingDataColor = [0.9 0.9 0.9];
    h1.MissingDataLabel = "n/a";
    h1.FontSize = 14;
    h1.Title = 'Outputs - R1/R2/N1';
    h1.XLabel = 'VEGFconc (ng/ml)';
    h1.YLabel = 'PlGFconc (ng/ml)';
    a = [0 ngml;ngml' a'];
    csvname = strcat("outputDataDr4/",calcname,"a.csv");
    csvwrite(csvname,a);
      ax2 = subplot(4,1,2);
      h2 = heatmap(ngml,ngml,b'); % (:,1,:) vegf and timepoints; % (1,:,:) plgf and timepoints
      h2.Colormap = map1;
      h2.ColorLimits = [0,scalemax];
      h2.CellLabelFormat = '%.0f';
      h2.MissingDataColor = [0.9 0.9 0.9];
      h2.MissingDataLabel = "n/a";
      h2.FontSize = 14;
      h2.Title = 'Outputs - --/R2/N1';
      h2.XLabel = 'VEGFconc (ng/ml)';
      h2.YLabel = 'PlGFconc (ng/ml)';
      b = [0 ngml;ngml' b'];
      csvname = strcat("outputDataDr4/",calcname,"b.csv");
      csvwrite(csvname,b);
        ax3 = subplot(4,1,3);
        h3 = heatmap(ngml,ngml,c'); % (:,1,:) vegf and timepoints; % (1,:,:) plgf and timepoints
        h3.Colormap = map1;
        h3.ColorLimits = [0,scalemax];
    	h3.CellLabelFormat = '%.0f';
        h3.MissingDataColor = [0.9 0.9 0.9];
        h3.MissingDataLabel = "n/a";
        h3.FontSize = 14;
        h3.Title = 'Outputs - R1/--/N1';
        h3.XLabel = 'VEGFconc (ng/ml)';
        h3.YLabel = 'PlGFconc (ng/ml)';
        c = [0 ngml;ngml' c'];
        csvname = strcat("outputDataDr4/",calcname,"c.csv");
        csvwrite(csvname,c);
          ax4 = subplot(4,1,4);
          h4 = heatmap(ngml,ngml,d'); % (:,1,:) vegf and timepoints; % (1,:,:) plgf and timepoints
          h4.Colormap = map1;
          h4.ColorLimits = [0,scalemax];
      	  h4.CellLabelFormat = '%.0f';
          h4.MissingDataColor = [0.9 0.9 0.9];
          h4.MissingDataLabel = "n/a";
          h4.FontSize = 14;
          h4.Title = 'Outputs - R1/R2/--';
          h4.XLabel = 'VEGFconc (ng/ml)';
          h4.YLabel = 'PlGFconc (ng/ml)';
          d = [0 ngml;ngml' d'];
          csvname = strcat("outputDataDr4/",calcname,"d.csv");
          csvwrite(csvname,d);
    set(eval(figname), 'Position',[0 0 1000 1000])
    set(eval(figname), 'PaperUnits', 'inches');
    set(eval(figname), 'PaperSize', [4 4]);
    exportgraphics(eval(figname),strcat("outputFigsDr4/",calcname,".png"),'Resolution',300);
end


%% Function - Run VEGF x PlGF simulation array

function [Routs] = RunVEGFPlGFarray(runModel,speciesInit,parameters,lig1,lig2,lig1conc,lig2conc,timepts)

if lig1 == "v165a"
    v165a = lig1conc;   
    v121a = v165a*0;
elseif lig1 == "v121a"
    v121a = lig1conc;
    v165a = v121a*0;
end

if lig2 == "plgf1"
    plgf1 = lig2conc;
    plgf2 = plgf1*0;
elseif lig2 == "plgf2"
    plgf2 = lig2conc;
    plgf1 = plgf2*0;
end

Routs.R1_total = zeros(length(lig1conc),length(lig2conc),length(timepts));
Routs.R2_total = Routs.R1_total;
Routs.N1_total = Routs.R1_total;
Routs.R1_surf = Routs.R1_total;
Routs.R2_surf = Routs.R1_total;
Routs.N1_surf = Routs.R1_total;
Routs.R1_int = Routs.R1_total;
Routs.R2_int = Routs.R1_total;
Routs.N1_int = Routs.R1_total;
Routs.R1_Rab4 = Routs.R1_total;
Routs.R2_Rab4 = Routs.R1_total;
Routs.N1_Rab4 = Routs.R1_total;
Routs.R1_Rab11 = Routs.R1_total;
Routs.R2_Rab11 = Routs.R1_total;
Routs.N1_Rab11 = Routs.R1_total;
Routs.R1R1_165a_total = Routs.R1_total;
Routs.R2R2_165a_total = Routs.R1_total;
Routs.R1R1_165a_surf = Routs.R1_total;
Routs.R2R2_165a_surf = Routs.R1_total;
Routs.R1R1_165a_int = Routs.R1_total;
Routs.R2R2_165a_int = Routs.R1_total;
Routs.R1R1_plgf1_total = Routs.R1_total;
Routs.R1R1_plgf1_surf = Routs.R1_total;
Routs.R1R1_plgf1_int = Routs.R1_total;

Routs.R1R1_165a_total_ACT = Routs.R1_total;
Routs.R2R2_165a_total_ACT = Routs.R1_total;
Routs.R1R1_165a_surf_ACT = Routs.R1_total;
Routs.R2R2_165a_surf_ACT = Routs.R1_total;
Routs.R1R1_165a_int_ACT = Routs.R1_total;
Routs.R2R2_165a_int_ACT = Routs.R1_total;
Routs.R1R1_plgf1_total_ACT = Routs.R1_total;
Routs.R1R1_plgf1_surf_ACT = Routs.R1_total;
Routs.R1R1_plgf1_int_ACT = Routs.R1_total;

Routs.R1R1_121a_total_ACT = Routs.R1_total;
Routs.R2R2_121a_total_ACT = Routs.R1_total;
Routs.R1R1_121a_surf_ACT = Routs.R1_total;
Routs.R2R2_121a_surf_ACT = Routs.R1_total;
Routs.R1R1_121a_int_ACT = Routs.R1_total;
Routs.R2R2_121a_int_ACT = Routs.R1_total;
Routs.R1R1_plgf2_total_ACT = Routs.R1_total;
Routs.R1R1_plgf2_surf_ACT = Routs.R1_total;
Routs.R1R1_plgf2_int_ACT = Routs.R1_total;

Routs.V165_surf = Routs.R1_total;
Routs.V165_int = Routs.R1_total;
Routs.V165_Rab4 = Routs.R1_total;
Routs.V165_Rab11 = Routs.R1_total;
Routs.V121_surf = Routs.R1_total;
Routs.V121_int = Routs.R1_total;
Routs.V121_Rab4 = Routs.R1_total;
Routs.V121_Rab11 = Routs.R1_total;
Routs.PGF1_surf = Routs.R1_total;
Routs.PGF1_int = Routs.R1_total;
Routs.PGF1_Rab4 = Routs.R1_total;
Routs.PGF1_Rab11 = Routs.R1_total;
Routs.PGF2_surf = Routs.R1_total;
Routs.PGF2_int = Routs.R1_total;
Routs.PGF2_Rab4 = Routs.R1_total;
Routs.PGF2_Rab11 = Routs.R1_total;
Routs.ligand_surf = Routs.R1_total;
Routs.ligand_int = Routs.R1_total;

Routs.R1_empty_surf = Routs.R1_total;
Routs.R2_empty_surf = Routs.R1_total;
Routs.N1_empty_surf = Routs.R1_total;
Routs.R1_empty_Rab4 = Routs.R1_total;
Routs.R2_empty_Rab4 = Routs.R1_total;
Routs.N1_empty_Rab4 = Routs.R1_total;
Routs.R1_empty_Rab11 = Routs.R1_total;
Routs.R2_empty_Rab11 = Routs.R1_total;
Routs.N1_empty_Rab11 = Routs.R1_total;

for i = 1:length(lig1conc)
    for j = 1:length(lig2conc)
        [Routs.R1_total(i,j,:), Routs.R2_total(i,j,:), Routs.N1_total(i,j,:), Routs.R1_surf(i,j,:), Routs.R2_surf(i,j,:), Routs.N1_surf(i,j,:), Routs.R1_int(i,j,:), Routs.R2_int(i,j,:), Routs.N1_int(i,j,:), Routs.R1_Rab4(i,j,:), Routs.R2_Rab4(i,j,:), Routs.N1_Rab4(i,j,:), Routs.R1_Rab11(i,j,:), Routs.R2_Rab11(i,j,:), Routs.N1_Rab11(i,j,:),...
            Routs.R1R1_165a_total(i,j,:), Routs.R2R2_165a_total(i,j,:), Routs.R1R1_165a_surf(i,j,:), Routs.R2R2_165a_surf(i,j,:), Routs.R1R1_165a_int(i,j,:), Routs.R2R2_165a_int(i,j,:), Routs.R1R1_plgf1_total(i,j,:), Routs.R1R1_plgf1_surf(i,j,:), Routs.R1R1_plgf1_int(i,j,:), ...
            Routs.R1R1_165a_total_ACT(i,j,:), Routs.R2R2_165a_total_ACT(i,j,:), Routs.R1R1_165a_surf_ACT(i,j,:), Routs.R2R2_165a_surf_ACT(i,j,:), Routs.R1R1_165a_int_ACT(i,j,:), Routs.R2R2_165a_int_ACT(i,j,:), Routs.R1R1_plgf1_total_ACT(i,j,:), Routs.R1R1_plgf1_surf_ACT(i,j,:), Routs.R1R1_plgf1_int_ACT(i,j,:), ...
            Routs.R1R1_121a_total_ACT(i,j,:), Routs.R2R2_121a_total_ACT(i,j,:), Routs.R1R1_121a_surf_ACT(i,j,:), Routs.R2R2_121a_surf_ACT(i,j,:), Routs.R1R1_121a_int_ACT(i,j,:), Routs.R2R2_121a_int_ACT(i,j,:), Routs.R1R1_plgf2_total_ACT(i,j,:), Routs.R1R1_plgf2_surf_ACT(i,j,:), Routs.R1R1_plgf2_int_ACT(i,j,:), ...
            Routs.V165_surf(i,j,:),Routs.V165_int(i,j,:),Routs.V165_Rab4(i,j,:),Routs.V165_Rab11(i,j,:),Routs.V121_surf(i,j,:),Routs.V121_int(i,j,:),Routs.V121_Rab4(i,j,:),Routs.V121_Rab11(i,j,:), ...
            Routs.PGF1_surf(i,j,:),Routs.PGF1_int(i,j,:),Routs.PGF1_Rab4(i,j,:),Routs.PGF1_Rab11(i,j,:),Routs.PGF2_surf(i,j,:),Routs.PGF2_int(i,j,:),Routs.PGF2_Rab4(i,j,:),Routs.PGF2_Rab11(i,j,:), Routs.ligand_surf(i,j,:), Routs.ligand_int(i,j,:),...
            Routs.R1_empty_surf(i,j,:),Routs.R2_empty_surf(i,j,:),Routs.N1_empty_surf(i,j,:),Routs.R1_empty_Rab4(i,j,:),Routs.R2_empty_Rab4(i,j,:),Routs.N1_empty_Rab4(i,j,:),Routs.R1_empty_Rab11(i,j,:),Routs.R2_empty_Rab11(i,j,:),Routs.N1_empty_Rab11(i,j,:)] ...
            = model20230914(runModel,speciesInit,parameters,v165a(i),v121a(i),plgf1(j),plgf2(j),timepts);
    end
    k = (i-1)*length(lig2conc)+j;
    disp(strcat('vegf x plgf',num2str(k*100/(length(lig1conc)*length(lig2conc))),'% complete'))
end

end


%% Function - Visualizing steady-state results for confirmation

function f1 = VizConcOut(ConcOut,timepoints)

ConcOutNames = fieldnames(ConcOut);
TotPanels = length(ConcOutNames);
PanelRows = ceil(length(ConcOutNames)/3);

f1 = figure('visible','off');
for m=1:TotPanels
	eval(append("ax",num2str(m),"=subplot(",num2str(PanelRows),",3,",num2str(m),");"));
	plot(timepoints/3600,ConcOut.(ConcOutNames{m}));
	hold on;
	title(ConcOutNames{m}, 'Interpreter','none');
end

end
