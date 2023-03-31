clear all;
close all;

% Driver for global and local sensitivity of unliganded trafficking parameters 
%  for VEGFR1, VEGFR2, NRP1
%  Rab7a is the degradation compartment

%% KEY INITIALIZATIONS and OPTIONS

model = "Unligated_VEGFR_model_20230321"; % which compiled model to use
baseparams = "SSfit2022.csv"; % which baseline unliganded receptor parameters to use
% Alternate parameter sets to use as a baseline: 'LWCfit2015.csv'
parameters = base_parameters(model,baseparams); % Initialize model parameters
% parameters.kN1R1on = 0; % if preventing R1-N1 coupling (e.g. Fig S5)
timestp_stst = 3600*3; % steady-state - only output every 3 hours
% time step size for the ligand simulations is defined in the ligand simulation code (RunOneVEGFPlGFsim.m) 

% Target values for surface receptor levels
R1surf_target = 1.8;
R2surf_target = 4.9;
N1surf_target = 68;
recsurftargets = [R1surf_target; R2surf_target; N1surf_target];

suppressligands = 1; % for unliganded simulations, set = 1 to suppress ligand outputs
% for liganded simulations, set = 0 to get both

RunSim = 2;
optimizeR = 0;

% Simulations: 
% 1 = simple global sensitivity - scan of each parameter (n points)
% 2 = simple local sensitivity - heatmap (many outputs, many inputs) - all run to st st
% 3 = simple local sensitivity - heatmap (many outputs, many inputs) - for 60 mins change


%% INITIAL BASELINE

% define global variables for tracking optimization loops (for debugging purposes)
global r1loop
global ploop
r1loop = 0;
ploop = 0;

%   Optimize production rates to surface receptor levels, if needed

if optimizeR == 1
    currentR1prodguess = parameters.kR1prod;
    currentR2prodguess = parameters.kR2prod;
    currentN1prodguess = parameters.kN1prod;
    p2(1)=currentR1prodguess;
    p2(2)=currentR2prodguess;
    p2(3)=currentN1prodguess; 
    initguesses = p2;
    lb=p2*.001;     % lower bound
    ub=p2*1000;     % upper bound
    options = optimoptions('lsqnonlin','Display','iter'); % display output
    optimalR = lsqnonlin(@CostFxn_ReceptorProduction,p2,lb,ub,options,recsurftargets,parameters,model);

    parameters.kR1prod = optimalR(1);
    parameters.kR2prod = optimalR(2);
    parameters.kN1prod = optimalR(3);
end

[timepoints, species_out, observables_out] = SimToSteadyState(parameters,model);
speciesInit = species_out(end,:); % we just ran everything to steady state, so use the outputs as new inputs
ConcOut=CalcOutputs(observables_out,suppressligands);
baselines = CalcBaselines(ConcOut);

ConcOutNames = fieldnames(ConcOut);
TotPanels = length(ConcOutNames);
PanelRows = ceil(length(ConcOutNames)/3);

Finalss = figure;
for m=1:TotPanels
    eval(append("ax",num2str(m),"=subplot(",num2str(PanelRows),",3,",num2str(m),");"));
    plot(timepoints/3600,ConcOut.(ConcOutNames{m}));
    hold on;
    title(ConcOutNames{m}, 'Interpreter','none');
end
        
%         createFigs = 0;
%         [R1_exp,R1_sim,R2_exp,R2_sim,N1_exp,N1_sim] = UnligandedHUVECExperimentsKQ(parameters,speciesInit,baselines,perturbs,model,createFigs);
        
p0.kR1prod      = parameters.kR1prod;
p0.kR1Rab5a     = parameters.kR1Rab5a; 
p0.kR1Rab4at7a  = parameters.kR1Rab4at7a;
p0.kR1Rab4a     = parameters.kR1Rab4a;
p0.kR1Rab4at11a = parameters.kR1Rab4at11a;
p0.kR1Rab11a    = parameters.kR1Rab11a;

p0.kR2prod      = parameters.kR2prod;
p0.kR2Rab5a     = parameters.kR2Rab5a;         
p0.kR2Rab4at7a  = parameters.kR2Rab4at7a;     
p0.kR2Rab4a     = parameters.kR2Rab4a;    
p0.kR2Rab4at11a = parameters.kR2Rab4at11a;     
p0.kR2Rab11a    = parameters.kR2Rab11a;

p0.kN1prod      = parameters.kN1prod;
p0.kN1Rab5a     = parameters.kN1Rab5a;
p0.kN1Rab4at7a  = parameters.kN1Rab4at7a;
p0.kN1Rab4a     = parameters.kN1Rab4a;
p0.kN1Rab4at11a = parameters.kN1Rab4at11a;
p0.kN1Rab11a    = parameters.kN1Rab11a;

p0names = fieldnames(p0);
expsimnames = {'R1CHX30m','R1CHX60m','R1CHX90m','R1CHX120m','R1CHX240m','R1siRab4','R1siRab11','R1siRab411','R1surf',...
               'R2CHX30m','R2CHX60m','R2CHX90m','R2CHX120m','R2CHX240m','R2siRab4','R2siRab11','R2siRab411','R2surf',...
               'N1CHX30m','N1CHX60m','N1CHX90m','N1CHX120m','N1CHX240m','N1siRab4','N1siRab11','N1siRab411','N1surf'}';

switch RunSim

%% GLOBAL SENSITIVITY CASES

    case 1

        scanparams = logspace(-2,2,17);
        for i = 1:length(p0names) % do all 15 parameters
            runparams = parameters;
            for j = 1:length(scanparams)
                ploop = i*100+j;
                r1loop = 0;
                runparams.(p0names{i}) = parameters.(p0names{i})*scanparams(j);
                runparams.kR1N1Rab5a     = runparams.kR1Rab5a; 
                runparams.kR1N1Rab4at7a  = runparams.kR1Rab4at7a;
                runparams.kR1N1Rab4a     = runparams.kR1Rab4a;
                runparams.kR1N1Rab4at11a = runparams.kR1Rab4at11a;
                runparams.kR1N1Rab11a    = runparams.kR1Rab11a;

%                [RUN SIMULATION pt 1 - steady state - GET DATA]
                if optimizeR == 1
                    currentR1prodguess = parameters.kR1prod;
                    currentR2prodguess = parameters.kR2prod;
                    currentN1prodguess = parameters.kN1prod;
                    p2(1)=currentR1prodguess;
                    p2(2)=currentR2prodguess;
                    p2(3)=currentN1prodguess; 
                    initguesses = p2;
                    lb=p2*.001;     % lower bound
                    ub=p2*1000;     % upper bound
                    options = optimoptions('lsqnonlin','Display','iter'); % display output
                    optimalR = lsqnonlin(@CostFxn_ReceptorProduction,p2,lb,ub,options,recsurftargets,runparams,model);

                    runparams.kR1prod = optimalR(1);
                    runparams.kR2prod = optimalR(2);
                    runparams.kN1prod = optimalR(3);
                end

                [timepoints, species_out, observables_out] = SimToSteadyState(runparams,model);
                speciesInit = species_out(end,:); % we just ran everything to steady state, so use the outputs as new inputs
                ConcOut=CalcOutputs(observables_out,suppressligands);
                baselines = CalcBaselines(ConcOut);
%                [DATA: surfperc x 3; st st R surf level x 3];
                surfPercR1.(p0names{i})(j) = ConcOut.R1_percSurf(end);
                surfPercR2.(p0names{i})(j) = ConcOut.R2_percSurf(end);
                surfPercN1.(p0names{i})(j) = ConcOut.N1_percSurf(end);
                surfR1.(p0names{i})(j) = ConcOut.R1_surf(end);
                surfR2.(p0names{i})(j) = ConcOut.R2_surf(end);
                surfN1.(p0names{i})(j) = ConcOut.N1_surf(end);
                totalR1.(p0names{i})(j) = ConcOut.R1_total(end);
                totalR2.(p0names{i})(j) = ConcOut.R2_total(end);
                totalN1.(p0names{i})(j) = ConcOut.N1_total(end);
                
%                [RUN SIMULATION pt 2 - KQ cost - GET DATA]
        % RUN KQ Unliganded EXPERIMENTS and output figures
        % severity of knockdowns (can put these in function call to see the effect of varying effect size)
                perturbs.R1reductionByCHX = 0.00; % 0 = complete shutdown of production
                perturbs.R2reductionByCHX = 0.00; % 0.5 = half of production
                perturbs.N1reductionByCHX = 0.00; % 1.0 = no change to production

                perturbs.Rab4reductionBysiRNA    = 1 - 0.80;  %80% reduction
                perturbs.Rab11reductionBysiRNA   = 1 - 0.80;

                createFigs = 0;
                [R1_exp,R1_sim,R2_exp,R2_sim,N1_exp,N1_sim] = UnligandedHUVECExperimentsAll(runparams,speciesInit,baselines,perturbs,recsurftargets,model,createFigs);
% Traditional CostOut - absolute differences in levels
%                 CostOut = [R1_sim - R1_exp; R2_sim - R2_exp; N1_sim - N1_exp]; %3*length(R1_theor) elements

% Consider alternate cost methods - (a) normalization
                CostOut = [(R1_sim - R1_exp)./R1_exp; (R2_sim - R2_exp)./R2_exp; (N1_sim - N1_exp)./N1_exp]; %3*length(R1_theor) elements
                CostOut (isnan(CostOut)) = 0; 

% Consider alternate cost methods - (b) fold change - max increase or decrease
%                 CostOuta = [R1_sim./R1_exp; R2_sim./R2_exp; N1_sim./N1_exp]; %3*length(R1_theor) elements
%                 CostOuta(isnan(CostOuta)) = 0;
%                 CostOutb = [R1_exp./R1_sim; R2_exp./R2_sim; N1_exp./N1_sim]; %3*length(R1_theor) elements
%                 CostOutb(isnan(CostOutb)) = 0;
%                 CostOut = max(CostOuta,CostOutb);

                k=9;
                CostOut(k) = CostOut(k)*10; % give these equal weight to the other readouts; basically, try to enforce low surface R1, 17 is the number of data points in Bouchert et al (2017) paper
                CostOut(k+length(R1_sim)) = CostOut(k+length(R1_sim))*10; % give these equal weight to the other readouts; basically, try to enforce low surface R1
                CostOut(k+length(R1_sim)+length(R2_sim)) = CostOut(k+length(R1_sim)+length(R2_sim))*10; % give these equal weight to the other readouts; basically, try to enforce low surface R1

                % Emphasize the siRNA simulations
                k=6:8;
                CostOut(k) = CostOut(k)*3; 
                CostOut(k+length(R1_sim)) = CostOut(k+length(R1_sim))*3; 
                CostOut(k+length(R1_sim)+length(R2_sim)) = CostOut(k+length(R1_sim)+length(R2_sim))*3; 
 
                costTotal.(p0names{i})(j) = sqrt(sum(CostOut.*CostOut));
                k = [9 9+length(R1_sim) 9+length(R1_sim)+length(R2_sim)];
                costPercSurf.(p0names{i})(j) = sqrt(sum(CostOut(k).*CostOut(k)));
                k = [1:5 1:5+length(R1_sim) 1:5+length(R1_sim)+length(R2_sim)];
                costCHX.(p0names{i})(j) = sqrt(sum(CostOut(k).*CostOut(k)));
                k = [6:8 6:8+length(R1_sim) 6:8+length(R1_sim)+length(R2_sim)];
                costRab411.(p0names{i})(j) = sqrt(sum(CostOut(k).*CostOut(k)));

            disp(['global univariate sensitivity ' num2str((length(scanparams)*(i-1)+j)/(length(scanparams)*length(p0names))*100) ' % complete']);

            end
            
     
        
            f1 = figure;
            ax1 = subplot(4,3,1);
                plot(scanparams,surfPercR1.(p0names{i})); set(gca,'Xscale','log');
                ylabel('R1 surface %'); xlabel('parameter value/baseline value'); 
            ax2 = subplot(4,3,2);
                plot(scanparams,surfPercR2.(p0names{i})); set(gca,'Xscale','log');
                ylabel('R2 surface %'); xlabel('parameter value/baseline value');
            ax3 = subplot(4,3,3);
                plot(scanparams,surfPercN1.(p0names{i})); set(gca,'Xscale','log');
                ylabel('N1 surface %'); xlabel('parameter value/baseline value');
            ax4 = subplot(4,3,4);
                plot(scanparams,surfR1.(p0names{i})); set(gca,'Xscale','log');
                ylabel('R1 surface #'); xlabel('parameter value/baseline value');
            ax5 = subplot(4,3,5);
                plot(scanparams,surfR2.(p0names{i})); set(gca,'Xscale','log');
                ylabel('R2 surface #'); xlabel('parameter value/baseline value');
            ax6 = subplot(4,3,6);
                plot(scanparams,surfN1.(p0names{i})); set(gca,'Xscale','log');
                ylabel('N1 surface #'); xlabel('parameter value/baseline value');
            ax7 = subplot(4,3,7);
                plot(scanparams,totalR1.(p0names{i})); set(gca,'Xscale','log');
                ylabel('R1 total #'); xlabel('parameter value/baseline value');
            ax8 = subplot(4,3,8);
                plot(scanparams,totalR2.(p0names{i})); set(gca,'Xscale','log');
                ylabel('R2 total #'); xlabel('parameter value/baseline value');
            ax9 = subplot(4,3,9);
                plot(scanparams,totalN1.(p0names{i})); set(gca,'Xscale','log');
                ylabel('N1 total #'); xlabel('parameter value/baseline value');
            ax10 = subplot(4,3,10);
                plot(scanparams,costTotal.(p0names{i})); set(gca,'Xscale','log'); set(gca,'Yscale','log');
                ylabel('costfxn total'); xlabel('parameter value/baseline value');
            ax11 = subplot(4,3,11);
                plot(scanparams,costCHX.(p0names{i})); set(gca,'Xscale','log'); set(gca,'Yscale','log');
                ylabel('costfxn chx'); xlabel('parameter value/baseline value');
            ax12 = subplot(4,3,12);
                plot(scanparams,costRab411.(p0names{i})); set(gca,'Xscale','log'); set(gca,'Yscale','log');
                ylabel('costfxn rab411'); xlabel('parameter value/baseline value');
            sgtitle (p0names{i});
%           add cost function results [
%           output plot results to files

        end
 
        a = zeros(length(p0names),length(scanparams));
        b = a; c = a; d = a; e = a;
        a2 = a; b2 = a; c2 = a; d2 = a; e2 = a;
        for i = 1:length(p0names)
        for j = 1:length(scanparams)
            a0(i,j) = costTotal.(p0names{i})(j);            
            b0(i,j) = costPercSurf.(p0names{i})(j);            
            c0(i,j) = costCHX.(p0names{i})(j);            
            d0(i,j) = costRab411.(p0names{i})(j);            
            e0(i,j) = surfPercR1.(p0names{i})(j);            
            f0(i,j) = surfPercR2.(p0names{i})(j);            
        end
        end
        a = a0/(min(min(a0))); b = b0/(min(min(b0)));
        c = c0/(min(min(c0))); d = d0/(min(min(d0))); 
        e = e0*100; f = f0*100;

        a2 = log10(a); b2 = log10(b);
        c2 = log10(c); d2 = log10(d); 
        e2 = log10(e); f2 = log10(f);        

        
        f2 = figure;
        for i = 1:length(p0names) % do all 15 parameters
            ax1 = subplot(3,6,i);
                plot(scanparams,costTotal.(p0names{i})); set(gca,'Xscale','log'); ylim([10^(floor(log10(min(min(a0))))) 10^(ceil(log10(max(max(a0)))))]); set(gca,'Yscale','log');
                ylabel('costfxn total'); xlabel('param(/baseline)');
                title(p0names{i});                
        end
        sgtitle ('cost fxn total');

        set(f2, 'Position',[0 0 1000 1000])
        set(f2, 'PaperUnits', 'inches');
        set(f2, 'PaperSize', [4 4]);
        exportgraphics(f2,strcat("GlobalSensitivity_CostTotal_RSS.png"),'Resolution',300);

        f3 = figure;
        for i = 1:length(p0names) % do all 15 parameters
            ax1 = subplot(3,6,i);
                plot(scanparams,surfPercR1.(p0names{i})); set(gca,'Xscale','log'); ylim([min(min(e0)) max(max(e0))]); % set(gca,'Yscale','log'); 
                ylabel('R1 surface %'); xlabel('param(/baseline)');
                title(p0names{i});                
        end
        sgtitle ('R1 surface percentage total');

        set(f3, 'Position',[0 0 1000 1000])
        set(f3, 'PaperUnits', 'inches');
        set(f3, 'PaperSize', [4 4]);
        exportgraphics(f3,strcat("GlobalSensitivity_R1surfperc.png"),'Resolution',300);

        f4 = figure;
        for i = 1:length(p0names) % do all 15 parameters
            ax1 = subplot(3,6,i);
                plot(scanparams,surfPercR2.(p0names{i})); set(gca,'Xscale','log'); ylim([min(min(f0)) max(max(f0))]); % set(gca,'Yscale','log'); 
                ylabel('R2 surface %'); xlabel('param(/baseline)');
                title(p0names{i});                
        end
        sgtitle ('R2 surface percentage total');

        set(f4, 'Position',[0 0 1000 1000])
        set(f4, 'PaperUnits', 'inches');
        set(f4, 'PaperSize', [4 4]);
        exportgraphics(f4,strcat("GlobalSensitivity_R2surfperc.png"),'Resolution',300);

        


        fs1 = figure;
        ax1 = subplot(4,1,1);
        h1 = htmp2(round(scanparams,2,'significant'),p0names,a2,'Cost total (log10 RSS)','parameter value/baseline value','Params',[min(min(a2)) max(max(a2))]);
            ax2 = subplot(4,1,2);
            h2 = htmp2(round(scanparams,2,'significant'),p0names,b2,'Cost PercSurf (log10 RSS)','parameter value/baseline value','Params',[min(min(b2)) max(max(b2))]);
                ax3 = subplot(4,1,3);
                h3 = htmp2(round(scanparams,2,'significant'),p0names,c2,'Cost CHX (log10 RSS)','parameter value/baseline value','Params',[min(min(c2)) max(max(c2))]);
                    ax4 = subplot(4,1,4);
                    h4 = htmp2(round(scanparams,2,'significant'),p0names,d2,'Cost Rab411 (log10 RSS)','parameter value/baseline value','Params',[min(min(d2)) max(max(d2))]);

        set(fs1, 'Position',[0 0 500 2000])
        set(fs1, 'PaperUnits', 'inches');
        set(fs1, 'PaperSize', [4 6]);
        exportgraphics(fs1,strcat("GlobalSensitivityHM1.png"),'Resolution',300);

        fs2 = figure;
        ax1 = subplot(4,1,1);
        h1 = htmp2(round(scanparams,2,'significant'),p0names,a,'Cost total (RSS)','parameter value/baseline value','Params',[min(min(a)) max(max(a))]);
            ax2 = subplot(4,1,2);
            h2 = htmp2(round(scanparams,2,'significant'),p0names,b,'Cost PercSurf (RSS)','parameter value/baseline value','Params',[min(min(b)) max(max(b))]);
                ax3 = subplot(4,1,3);
                h3 = htmp2(round(scanparams,2,'significant'),p0names,c,'Cost CHX (RSS)','parameter value/baseline value','Params',[min(min(c)) max(max(c))]);
                    ax4 = subplot(4,1,4);
                    h4 = htmp2(round(scanparams,2,'significant'),p0names,d,'Cost Rab411 (RSS)','parameter value/baseline value','Params',[min(min(d)) max(max(d))]);

        set(fs2, 'Position',[0 0 500 2000])
        set(fs2, 'PaperUnits', 'inches');
        set(fs2, 'PaperSize', [4 6]);
        exportgraphics(fs2,strcat("GlobalSensitivityHM2.png"),'Resolution',300);

        fs3 = figure;
        ax1 = subplot(2,1,1);
        h1 = htmp2(round(scanparams,2,'significant'),p0names,e,'PercSurfR1 (%)','parameter value/baseline value','Params',[min(min(e)) max(max(e))]);
            ax2 = subplot(2,1,2);
            h2 = htmp2(round(scanparams,2,'significant'),p0names,f,'PercSurfR2 (%)','parameter value/baseline value','Params',[min(min(f)) max(max(f))]);

        set(fs3, 'Position',[0 0 500 2000])
        set(fs3, 'PaperUnits', 'inches');
        set(fs3, 'PaperSize', [4 6]);
        exportgraphics(fs3,strcat("GlobalSensitivityHM3.png"),'Resolution',300);

        fs4 = figure;
        ax1 = subplot(2,1,1);
        h1 = htmp2(round(scanparams,2,'significant'),p0names,a,'Cost total (RSS)','parameter value/baseline value','Params',[min(min(a)) max(max(a))]);
            ax1 = subplot(2,1,2);
            h2 = htmp2(round(scanparams,2,'significant'),p0names,a2,'Cost total (log10 RSS)','parameter value/baseline value','Params',[min(min(a2)) max(max(a2))]);

        set(fs4, 'Position',[0 0 500 2000])
        set(fs4, 'PaperUnits', 'inches');
        set(fs4, 'PaperSize', [4 6]);
        exportgraphics(fs4,strcat("GlobalSensitivityHM4.png"),'Resolution',300);

        reframe = [2 8 14 3 9 15 4 10 16 5 11 17 6 12 18 1 7 13];
        a2reframed = a2(reframe,:);
        p0namesreframed = p0names(reframe);

        fs5 = figure;
        ax1 = subplot(2,1,1);
        h1 = htmp2(round(scanparams,2,'significant'),p0names,a2,'Cost total (log10 RSS)','parameter value/baseline value','Params',[min(min(a2)) max(max(a2))]);
            ax2 = subplot(2,1,2);
            h2 = htmp2(round(scanparams,2,'significant'),p0namesreframed,a2reframed,'Cost total (log10 RSS)','parameter value/baseline value','Params',[min(min(a2)) max(max(a2))]);

        set(fs5, 'Position',[0 0 500 800])
        set(fs5, 'PaperUnits', 'inches');
        set(fs5, 'PaperSize', [4 6]);
        exportgraphics(fs5,strcat("GlobalSensitivityHM5.png"),'Resolution',300);


% plot data
% loop through the parameters (3x5 or 3x6) and give total cost with ylog
% label

%% LOCAL SENSITIVITY CASES

    case 2

        ConcOut_bs=ConcOut;        
        ConcOut_bs_red = ConcOutReduced(ConcOut_bs);        
        ConcOut_bs_red2 = ConcOutReduced2(ConcOut_bs);        
        ConcOutNames = fieldnames(ConcOut_bs);        
        ConcOutRedNames = fieldnames(ConcOut_bs_red);        
        ConcOutRedNames2 = fieldnames(ConcOut_bs_red2);        
 
        % based on steady state, now run experiments/perturbations
        % RUN Unliganded EXPERIMENTS, quantify sensitivity of outputs
        % also need baseline simulation!
        % severity of knockdowns (can put these in function call to see the effect of varying effect size)
        perturbs.R1reductionByCHX = 0.00; % 0 = complete shutdown of production
        perturbs.R2reductionByCHX = 0.00; % 0.5 = half of production
        perturbs.N1reductionByCHX = 0.00; % 1.0 = no change to production

        perturbs.Rab4reductionBysiRNA    = 1 - 0.80;  %80% reduction
        perturbs.Rab11reductionBysiRNA   = 1 - 0.80;
        % set up outputs to show (calc outputs) - experimental data

        createFigs = 0;
        [R1_exp,R1_sim,R2_exp,R2_sim,N1_exp,N1_sim] = UnligandedHUVECExperimentsAll(parameters,speciesInit,baselines,perturbs,recsurftargets,model,createFigs);
        outputs_bs = [R1_sim(1:8)/baselines.R1_total;R1_sim(9)/baselines.R1_surf;...
                        R2_sim(1:8)/baselines.R2_total;R2_sim(9)/baselines.R2_surf;...
                        N1_sim(1:8)/baselines.N1_total;N1_sim(9)/baselines.N1_surf];

        ParamDelta = 1.05;
        
        results = zeros (length(p0names),length(ConcOutNames));
        sens = results;
        results_red = zeros (length(p0names),length(ConcOutRedNames));
        sens_red = results_red;
        results_red2 = zeros (length(p0names),length(ConcOutRedNames2));
        sens_red2 = results_red2;
        for i = 1:length(p0names)
            runparams = parameters;
            runparams.(p0names{i}) = parameters.(p0names{i})*ParamDelta;
                runparams.kR1N1Rab5a     = runparams.kR1Rab5a; 
                runparams.kR1N1Rab4at7a  = runparams.kR1Rab4at7a;
                runparams.kR1N1Rab4a     = runparams.kR1Rab4a;
                runparams.kR1N1Rab4at11a = runparams.kR1Rab4at11a;
                runparams.kR1N1Rab11a    = runparams.kR1Rab11a;

                 if optimizeR == 1
                   	ploop = i;
                    r1loop = 0;
                    currentR1prodguess = parameters.kR1prod;
                    currentR2prodguess = parameters.kR2prod;
                    currentN1prodguess = parameters.kN1prod;
                    p2(1)=currentR1prodguess;
                    p2(2)=currentR2prodguess;
                    p2(3)=currentN1prodguess; 
                    initguesses = p2;
                    lb=p2*.001;     % lower bound
                    ub=p2*1000;     % upper bound
                    options = optimoptions('lsqnonlin','Display','iter'); % display output
                    optimalR = lsqnonlin(@CostFxn_ReceptorProduction,p2,lb,ub,options,recsurftargets,runparams,model);

                    runparams.kR1prod = optimalR(1);
                    runparams.kR2prod = optimalR(2);
                    runparams.kN1prod = optimalR(3);
                end

                [timepoints, species_out, observables_out] = SimToSteadyState(runparams,model);
                speciesInit = species_out(end,:); % we just ran everything to steady state, so use the outputs as new inputs
                ConcOut=CalcOutputs(observables_out,suppressligands);
                tempbaselines = CalcBaselines(ConcOut);
%                [DATA: surfperc x 3; st st R surf level x 3];
                for j = 1:length(ConcOutNames)
                    results(i,j) = ConcOut.(ConcOutNames{j})(end);
                    bsline = ConcOut_bs.(ConcOutNames{j})(end);
                    sens(i,j) = ((results(i,j) - bsline) / bsline) / (ParamDelta-1);
                end
                ConcOutRed = ConcOutReduced(ConcOut);
                for j = 1:length(ConcOutRedNames)
                    results_red(i,j) = ConcOutRed.(ConcOutRedNames{j})(end);
                    bsline = ConcOut_bs_red.(ConcOutRedNames{j})(end);
                    sens_red(i,j) = ((results_red(i,j) - bsline) / bsline) / (ParamDelta-1);
                end
                ConcOutRed2 = ConcOutReduced2(ConcOut);
                for j = 1:length(ConcOutRedNames2)
                    results_red2(i,j) = ConcOutRed2.(ConcOutRedNames2{j})(end);
                    bsline = ConcOut_bs_red2.(ConcOutRedNames2{j})(end);
                    sens_red2(i,j) = ((results_red2(i,j) - bsline) / bsline) / (ParamDelta-1);
                end
                disp(['local univariate sensitivity ' num2str(i/length(p0names)*100) ' % complete']);  
                
                % based on steady state, now run experiments/perturbations
                % RUN Unliganded EXPERIMENTS, quantify sensitivity of outputs
                % also need baseline simulation!
                % severity of knockdowns (can put these in function call to see the effect of varying effect size)
                perturbs.R1reductionByCHX = 0.00; % 0 = complete shutdown of production
                perturbs.R2reductionByCHX = 0.00; % 0.5 = half of production
                perturbs.N1reductionByCHX = 0.00; % 1.0 = no change to production

                perturbs.Rab4reductionBysiRNA    = 1 - 0.80;  %80% reduction
                perturbs.Rab11reductionBysiRNA   = 1 - 0.80;
                % set up outputs to show (calc outputs) - experimental data

                createFigs = 0;
                [R1_exp,R1_sim,R2_exp,R2_sim,N1_exp,N1_sim] = UnligandedHUVECExperimentsAll(runparams,speciesInit,tempbaselines,perturbs,recsurftargets,model,createFigs);
                outputs = [R1_sim;R2_sim;N1_sim];
                outputs = [R1_sim(1:8)/tempbaselines.R1_total;R1_sim(9)/tempbaselines.R1_surf;...
                        R2_sim(1:8)/tempbaselines.R2_total;R2_sim(9)/tempbaselines.R2_surf;...
                        N1_sim(1:8)/tempbaselines.N1_total;N1_sim(9)/tempbaselines.N1_surf];
                for j = 1:length(outputs)
                    exptlresults(i,j) = outputs(j);
                    bsline = outputs_bs(j);
                    exptlsens(i,j) = ((exptlresults(i,j) - bsline) / bsline) / (ParamDelta-1);
                end

        end
        sens_red2(abs(sens_red2)<1e-4)= 0;
        sens_red (abs(sens_red)<1e-4) = 0;
        sens     (abs(sens)<1e-4)     = 0;
        exptlsens     (abs(exptlsens)<1e-4)     = 0;

        % set up outputs to show (calc outputs) - steady state
        fs1 = figure;

        ax1 = subplot(2,1,1);
        h1 = htmp(ConcOutNames,p0names,results,'Results','Output','Params');
            ax2 = subplot(2,1,2);
            h2 = htmp11(ConcOutNames,p0names,sens,'Sensitivity','Output','Params');

        set(fs1, 'Position',[0 0 1000 1000])
        set(fs1, 'PaperUnits', 'inches');
        set(fs1, 'PaperSize', [4 4]);
        exportgraphics(fs1,strcat("Sensitivity1.png"),'Resolution',300);

        fs2 = figure;

        ax1 = subplot(2,1,1);
        h1 = htmp(ConcOutRedNames,p0names,results_red,'Results','Output','Params');
            ax2 = subplot(2,1,2);
            h2 = htmp11(ConcOutRedNames,p0names,sens_red,'Sensitivity','Output','Params');

        set(fs2, 'Position',[0 0 1000 1000])
        set(fs2, 'PaperUnits', 'inches');
        set(fs2, 'PaperSize', [4 4]);
        exportgraphics(fs2,strcat("Sensitivity2.png"),'Resolution',300);

        fs3 = figure;

        ax1 = subplot(2,1,1);
        h1 = htmp(ConcOutRedNames2,p0names,results_red2,'Results','Output','Params');
            ax2 = subplot(2,1,2);
            h2 = htmp11(ConcOutRedNames2,p0names,sens_red2,'Sensitivity','Output','Params');

        set(fs3, 'Position',[0 0 1000 1000])
        set(fs3, 'PaperUnits', 'inches');
        set(fs3, 'PaperSize', [4 4]);
        exportgraphics(fs3,strcat("Sensitivity3.png"),'Resolution',300);

        
        fs4 = figure;

        ax1 = subplot(2,1,1);
        h1 = htmp(expsimnames,p0names,exptlresults,'Results','Output','Params');
            ax2 = subplot(2,1,2);
            h2 = htmp11(expsimnames,p0names,exptlsens,'Sensitivity','Output','Params');

        set(fs4, 'Position',[0 0 1000 1000])
        set(fs4, 'PaperUnits', 'inches');
        set(fs4, 'PaperSize', [4 4]);
        exportgraphics(fs4,strcat("Sensitivity4.png"),'Resolution',300);
        
        csvwrite('Local_Sensitivity_12outputs.csv',sens_red);
            
        % based on steady state, now run experiments/perturbations
         % RUN KQ Unliganded EXPERIMENTS and output figures
        % severity of knockdowns (can put these in function call to see the effect of varying effect size)
        perturbs.R1reductionByCHX = 0.00; % 0 = complete shutdown of production
        perturbs.R2reductionByCHX = 0.00; % 0.5 = half of production
        perturbs.N1reductionByCHX = 0.00; % 1.0 = no change to production

        perturbs.Rab4reductionBysiRNA    = 1 - 0.80;  %80% reduction
        perturbs.Rab11reductionBysiRNA   = 1 - 0.80;
        % set up outputs to show (calc outputs) - experimental data
        

    case 3

        ConcOut_bs=ConcOut;        
        ConcOut_bs_red = ConcOutReduced(ConcOut_bs);        
        ConcOut_bs_red2 = ConcOutReduced2(ConcOut_bs);        
        ConcOutNames = fieldnames(ConcOut_bs);        
        ConcOutRedNames = fieldnames(ConcOut_bs_red);        
        ConcOutRedNames2 = fieldnames(ConcOut_bs_red2);        
        ParamDelta = 1.05;
        
        results = zeros (length(p0names),length(ConcOutNames));
        sens = results;
        results_red = zeros (length(p0names),length(ConcOutRedNames));
        sens_red = results_red;
        results_red2 = zeros (length(p0names),length(ConcOutRedNames2));
        sens_red2 = results_red2;
        for i = 1:length(p0names)
            runparams = parameters;
            runparams.(p0names{i}) = parameters.(p0names{i})*ParamDelta;
                runparams.kR1N1Rab5a     = runparams.kR1Rab5a; 
                runparams.kR1N1Rab4at7a  = runparams.kR1Rab4at7a;
                runparams.kR1N1Rab4a     = runparams.kR1Rab4a;
                runparams.kR1N1Rab4at11a = runparams.kR1Rab4at11a;
                runparams.kR1N1Rab11a    = runparams.kR1Rab11a;

                 if optimizeR == 1
                   	ploop = i;
                    r1loop = 0;
                    currentR1prodguess = parameters.kR1prod;
                    currentR2prodguess = parameters.kR2prod;
                    currentN1prodguess = parameters.kN1prod;
                    p2(1)=currentR1prodguess;
                    p2(2)=currentR2prodguess;
                    p2(3)=currentN1prodguess; 
                    initguesses = p2;
                    lb=p2*.001;     % lower bound
                    ub=p2*1000;     % upper bound
                    options = optimoptions('lsqnonlin','Display','iter'); % display output
                    optimalR = lsqnonlin(@CostFxn_ReceptorProduction,p2,lb,ub,options,recsurftargets,runparams,model);

                    runparams.kR1prod = optimalR(1);
                    runparams.kR2prod = optimalR(2);
                    runparams.kN1prod = optimalR(3);
                 end
                 
                 Timepoints  = 0:60:3600;
%                 [timepoints, species_out, observables_out] = SimToSteadyState(runparams,model,timestp_stst);
                [~, timepoints, species_out, observables_out] = eval(append(model,"(Timepoints', speciesInit, runparams, 1)"));
%                 speciesInit = species_out(end,:); % we just ran everything to steady state, so use the outputs as new inputs
                ConcOut=CalcOutputs(observables_out,suppressligands);
%                [DATA: surfperc x 3; st st R surf level x 3];
                for j = 1:length(ConcOutNames)
                    results(i,j) = ConcOut.(ConcOutNames{j})(end);
                    bsline = ConcOut_bs.(ConcOutNames{j})(end);
                    sens(i,j) = ((results(i,j) - bsline) / bsline) / (ParamDelta-1);
                end
                ConcOutRed = ConcOutReduced(ConcOut);
                for j = 1:length(ConcOutRedNames)
                    results_red(i,j) = ConcOutRed.(ConcOutRedNames{j})(end);
                    bsline = ConcOut_bs_red.(ConcOutRedNames{j})(end);
                    sens_red(i,j) = ((results_red(i,j) - bsline) / bsline) / (ParamDelta-1);
                end
                ConcOutRed2 = ConcOutReduced2(ConcOut);
                for j = 1:length(ConcOutRedNames2)
                    results_red2(i,j) = ConcOutRed2.(ConcOutRedNames2{j})(end);
                    bsline = ConcOut_bs_red2.(ConcOutRedNames2{j})(end);
                    sens_red2(i,j) = ((results_red2(i,j) - bsline) / bsline) / (ParamDelta-1);
                end
                disp(['local univariate sensitivity ' num2str(i/length(p0names)*100) ' % complete']);        
        end
        sens_red2(abs(sens_red2)<1e-4)= 0;
        sens_red (abs(sens_red)<1e-4) = 0;
        sens     (abs(sens)<1e-4)     = 0;
        % set up outputs to show (calc outputs) - steady state
        fs1 = figure;

        ax1 = subplot(2,1,1);
        h1 = htmp(ConcOutNames,p0names,results,'Results','Output','Params');
            ax2 = subplot(2,1,2);
            h2 = htmp11(ConcOutNames,p0names,sens,'Sensitivity','Output','Params');

        set(fs1, 'Position',[0 0 1000 1000])
        set(fs1, 'PaperUnits', 'inches');
        set(fs1, 'PaperSize', [4 4]);
        exportgraphics(fs1,strcat("Sensitivity1.png"),'Resolution',300);

        fs2 = figure;

        ax1 = subplot(2,1,1);
        h1 = htmp(ConcOutRedNames,p0names,results_red,'Results','Output','Params');
            ax2 = subplot(2,1,2);
            h2 = htmp11(ConcOutRedNames,p0names,sens_red,'Sensitivity','Output','Params');

        set(fs2, 'Position',[0 0 1000 1000])
        set(fs2, 'PaperUnits', 'inches');
        set(fs2, 'PaperSize', [4 4]);
        exportgraphics(fs2,strcat("Sensitivity2.png"),'Resolution',300);

        fs3 = figure;

        ax1 = subplot(2,1,1);
        h1 = htmp(ConcOutRedNames2,p0names,results_red2,'Results','Output','Params');
            ax2 = subplot(2,1,2);
            h2 = htmp11(ConcOutRedNames2,p0names,sens_red2,'Sensitivity','Output','Params');

        set(fs3, 'Position',[0 0 1000 1000])
        set(fs3, 'PaperUnits', 'inches');
        set(fs3, 'PaperSize', [4 4]);
        exportgraphics(fs3,strcat("Sensitivity3.png"),'Resolution',300);

        


end



%% USEFUL FUNCTIONS

function a = ConcOutReduced(b)

a.R1total=b.R1_total;        
a.R1surf=b.R1_surf;        
a.R1int=b.R1_int;        
a.R1percSurf=b.R1_percSurf;        
a.R2total=b.R2_total;        
a.R2surf=b.R2_surf;        
a.R2int=b.R2_int;        
a.R2percSurf=b.R2_percSurf;        
a.N1total=b.N1_total;        
a.N1surf=b.N1_surf;        
a.N1int=b.N1_int;        
a.N1percSurf=b.N1_percSurf;        

end

function a = ConcOutReduced2(b)

a.R1total=b.R1_total;        
a.R2total=b.R2_total;        
a.N1total=b.N1_total;        
a.R1surf=b.R1_surf;        
a.R2surf=b.R2_surf;        
a.N1surf=b.N1_surf;        
a.R1int=b.R1_int;        
a.R2int=b.R2_int;        
a.N1int=b.N1_int;        
a.R1percSurf=b.R1_percSurf;        
a.R2percSurf=b.R2_percSurf;        
a.N1percSurf=b.N1_percSurf;        

end
