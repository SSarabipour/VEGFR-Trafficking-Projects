clear all;
close all;

% Driver for single-runs for given simulation parameters 
% VEGFR1, VEGFR2, NRP1 Trafficking Model 2022-2023

%% KEY INITIALIZATIONS and OPTIONS

model = "Unligated_VEGFR_model_20230831"; % which compiled model to use
baseparams = "SSfit2023.csv"; % which baseline unliganded receptor parameters to use
parameters = base_parameters(model,baseparams); % Initialize model parameters
% parameters.kN1R1on_surf = 0; % if preventing R1-N1 coupling (e.g. Fig S2)
% parameters.kN1R1on_rab4 = 0; % if preventing R1-N1 coupling (e.g. Fig S2)
% parameters.kN1R1on_rab11 = 0; % if preventing R1-N1 coupling (e.g. Fig S2)
% parameters.kN1prod = parameters.kN1prod/1000; % if preventing R1-N1 coupling (e.g. Fig S5)

timestp_stst = 3600*3; % steady-state - only output every 3 hours
% time step size for the ligand simulations is defined in the ligand simulation code (RunOneVEGFPlGFsim.m) 

% Target values for surface receptor levels
R1surf_target = 1800;
R2surf_target = 4900;
N1surf_target = 68000;
recsurftargets = [R1surf_target; R2surf_target; N1surf_target];

suppressligands = 1; % for unliganded simulations, set = 1 to suppress ligand outputs
% for liganded simulations, set = 0 to get both

RunSim = 1;
optimizeR = 0;
% Simulations: 
% 1 = single run to steady state PLUS unliganded trafficking experiments for a given set of parameters
% 2 = multiple runs to steady state PLUS unliganded trafficking experiments
% for a given set of parameters - to calculate costs
% 3 = single run to steady state PLUS CHQ simulations
% 4 = set of runs for different R1-N1 coupling rates
% 5 = set of runs for different R1-R1 coupling rates, no kcR1N1 if needed

switch RunSim

%% SINGLE RUNS

%% RUN UNLIGANDED TRAFFICKING PARAMETERS (Not optimizing) - for a specific set of trafficking parameters

    case 1

        % define global variables for tracking optimization loops (for debugging purposes)
        global r1loop
        global ploop
        r1loop = 0;
        ploop = 0;

% SELECT WHICH ONE PARAMETER SET TO RUN:        
%         p1 = readmatrix('LWCfit2015.csv');
%         p1 = readmatrix('SSfit2022.csv');
        p1 = readmatrix('SSfit2023.csv');
        
        parameters.kR1Rab5a     = p1(1); 
        parameters.kR1Rab4at7a  = p1(2); 
        parameters.kR1Rab4a     = p1(3); 
        parameters.kR1Rab4at11a = p1(4); 
        parameters.kR1Rab11a    = p1(5); 

        parameters.kR1N1Rab5a     = p1(1);
        parameters.kR1N1Rab4at7a  = p1(2); 
        parameters.kR1N1Rab4a     = p1(3); 
        parameters.kR1N1Rab4at11a = p1(4); 
        parameters.kR1N1Rab11a    = p1(5); 
        
        parameters.kR2Rab5a     = p1(6); 
        parameters.kR2Rab4at7a  = p1(7); 
        parameters.kR2Rab4a     = p1(8); 
        parameters.kR2Rab4at11a = p1(9); 
        parameters.kR2Rab11a    = p1(10);
        
        parameters.kN1Rab5a     = p1(11);
        parameters.kN1Rab4at7a  = p1(12);
        parameters.kN1Rab4a     = p1(13);
        parameters.kN1Rab4at11a = p1(14);
        parameters.kN1Rab11a    = p1(15);

        parameters.kR1prod = p1(16); 
        parameters.kR2prod = p1(17); 
        parameters.kN1prod = p1(18); 

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

% Alternative - particleswarm optimization (inefficient for this monotonic problem)
%             global PStargets
%             global PSparameters
%             global PSmodel
%             PStargets = recsurftargets;
%             PSparameters = parameters;
%             PSmodel = model;
%             options = optimoptions('particleswarm','SwarmSize',100,'HybridFcn',@fmincon,'Display','iter');
%             options = optimoptions('particleswarm','SwarmSize',100,'Display','iter');
%             [optimalR, fval, exitflag, output] = particleswarm(@CostFxn_ReceptorProductionPS,length(p2),lb,ub,options);

            parameters.kR1prod = optimalR(1);
            parameters.kR2prod = optimalR(2);
            parameters.kN1prod = optimalR(3);
        end

    % RUN Unliganded EXPERIMENTS and output figures
        % severity of knockdowns 
        perturbs.R1reductionByCHX = 0.00; % 0 = complete shutdown of production
        perturbs.R2reductionByCHX = 0.00; % 0.5 = half of production
        perturbs.N1reductionByCHX = 0.00; % 1.0 = no change to production

        perturbs.Rab4reductionBysiRNA    = 1 - 0.80;  %80% reduction
        perturbs.Rab11reductionBysiRNA   = 1 - 0.80;

        [timepoints, species_out, observables_out] = SimToSteadyState(parameters,model);
        speciesInit = species_out(end,:); % we just ran everything to steady state, so use the outputs as new inputs
        ConcOut=CalcOutputs(observables_out,suppressligands);
        baselines = CalcBaselines(ConcOut);
        perturbsnone = perturbs;
        perturbfields = fieldnames(perturbs);
        for k = 1:length(perturbfields)
            perturbsnone.(perturbfields{k}) = 1;
        end
        [R1rates,R2rates,N1rates] = CalcRates(ConcOut,parameters,perturbsnone);
        ratenames = fieldnames(R1rates);
        for k = 1:length(ratenames)
            R1ratesout(k) = R1rates.(ratenames{k})(end);
            R2ratesout(k) = R2rates.(ratenames{k})(end);
            N1ratesout(k) = N1rates.(ratenames{k})(end);
        end
        ratesout = [R1ratesout' R2ratesout' N1ratesout'];
        csvwrite('EndpointRates.csv',ratesout);
        concnames = fieldnames(baselines);
        for k = 1:length(concnames)
            concsout(k) = baselines.(concnames{k});
        end
        csvwrite('EndpointLevels.csv',concsout');
        
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
        
        createFigs = 1;
        [R1_exp,R1_sim,R2_exp,R2_sim,N1_exp,N1_sim] = UnligandedHUVECExperimentsAll(parameters,speciesInit,baselines,perturbs,recsurftargets,model,createFigs);

        
        
        
        
        
    case 2

        fits = readmatrix('100fits20230831.csv');
        inits = readmatrix('100fits20230831_inits.csv');
        fits = fits';
        inits = inits';
        FitRunCost = zeros(size(fits,2),1);
        FitR1_percSurf = FitRunCost;
        FitR2_percSurf = FitRunCost;
        FitN1_percSurf = FitRunCost;
        FitR1_Surf = FitRunCost;
        FitR2_Surf = FitRunCost;
        FitN1_Surf = FitRunCost;
        FitR1_Rab4 = FitRunCost;
        FitR2_Rab4 = FitRunCost;
        FitN1_Rab4 = FitRunCost;
        FitR1_Rab11 = FitRunCost;
        FitR2_Rab11 = FitRunCost;
        FitN1_Rab11 = FitRunCost;
        FitR1_Int = FitRunCost;
        FitR2_Int = FitRunCost;
        FitN1_Int = FitRunCost;
        FitR1_2hrCHX = FitRunCost;
        FitR2_2hrCHX = FitRunCost;
        FitN1_2hrCHX = FitRunCost;
        FitR1_siRNA = FitRunCost;
        FitR2_siRNA = FitRunCost;
        FitN1_siRNA = FitRunCost;

        for i = 1:size(fits,2) 
            p1 = fits(:,i);
                parameters.kR1Rab5a     = p1(1);
                parameters.kR1Rab4at7a  = p1(2);
                parameters.kR1Rab4a     = p1(3);
                parameters.kR1Rab4at11a = p1(4);
                parameters.kR1Rab11a    = p1(5);

                parameters.kR1N1Rab5a     = p1(1);
                parameters.kR1N1Rab4at7a  = p1(2); 
                parameters.kR1N1Rab4a     = p1(3); 
                parameters.kR1N1Rab4at11a = p1(4); 
                parameters.kR1N1Rab11a    = p1(5); 

                parameters.kR2Rab5a     = p1(6);         
                parameters.kR2Rab4at7a  = p1(7);
                parameters.kR2Rab4a     = p1(8);  
                parameters.kR2Rab4at11a = p1(9);     
                parameters.kR2Rab11a    = p1(10);   

                parameters.kN1Rab5a     = p1(11);
                parameters.kN1Rab4at7a  = p1(12);
                parameters.kN1Rab4a     = p1(13);
                parameters.kN1Rab4at11a = p1(14);   
                parameters.kN1Rab11a    = p1(15);

                parameters.kR1prod = p1(16);
                parameters.kR2prod = p1(17);
                parameters.kN1prod = p1(18);
            
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

            % RUN KQ Unliganded EXPERIMENTS and output figures
            % severity of knockdowns (can put these in function call to see the effect of varying effect size)
            perturbs.R1reductionByCHX = 0.00; % 0 = complete shutdown of production
            perturbs.R2reductionByCHX = 0.00; % 0.5 = half of production
            perturbs.N1reductionByCHX = 0.00; % 1.0 = no change to production

            perturbs.Rab4reductionBysiRNA    = 1 - 0.80;  %80% reduction
            perturbs.Rab11reductionBysiRNA   = 1 - 0.80;

            [timepoints, species_out, observables_out] = SimToSteadyState(parameters,model);
            speciesInit = species_out(end,:); % we just ran everything to steady state, so use the outputs as new inputs
            ConcOut=CalcOutputs(observables_out,suppressligands);
            baselines = CalcBaselines(ConcOut);
        
            ConcOutNames = fieldnames(ConcOut);
            TotPanels = length(ConcOutNames);
            PanelRows = length(ConcOutNames)/3;

% Visualize individual run results
%             Finalss = figure;
%             for m=1:TotPanels
%                 eval(append("ax",num2str(m),"=subplot(",num2str(PanelRows),",3,",num2str(m),");"));
%                 plot(timepoints/3600,ConcOut.(ConcOutNames{m}));
%                 hold on;
%                 title(ConcOutNames{m}, 'Interpreter','none');
%             end
            
            createFigs = 0;
            [R1_exp,R1_sim,R2_exp,R2_sim,N1_exp,N1_sim] = UnligandedHUVECExperimentsAll(parameters,speciesInit,baselines,perturbs,recsurftargets,model,createFigs);
            
            CostOut = [(R1_sim - R1_exp)./R1_exp; (R2_sim - R2_exp)./R2_exp; (N1_sim - N1_exp)./N1_exp]; %3*length(R1_theor) elements
            CostOut (isnan(CostOut)) = 0; 

            %% COST ADJUSTMENTS
            % For the surface-receptor-ratio specifically, we give the error more weight to help force compliance

%             % Emphasize the biotinylayion surface-internal simulations
%             j=9;
%             CostOut(j) = CostOut(j)*10; % give these equal weight to the other readouts; basically, try to enforce low surface R1, 17 is the number of data points in Bouchert et al (2017) paper
%             CostOut(j+length(R1_sim)) = CostOut(j+length(R1_sim))*10; % give these equal weight to the other readouts; basically, try to enforce low surface R1
%             CostOut(j+length(R1_sim)+length(R2_sim)) = CostOut(j+length(R1_sim)+length(R2_sim))*10; % give these equal weight to the other readouts; basically, try to enforce low surface R1
% 
%             % Emphasize the siRNA simulations
%             j=6:8;
%             CostOut(j) = CostOut(j)*3; 
%             CostOut(j+length(R1_sim)) = CostOut(j+length(R1_sim))*3; 
%             CostOut(j+length(R1_sim)+length(R2_sim)) = CostOut(j+length(R1_sim)+length(R2_sim))*3; 

            
%             FitRunCost(i) = sqrt(sum(CostOut.*CostOut));
            FitRunCost(i) = (sum(CostOut.*CostOut));
            FitR1_percSurf(i) = baselines.R1_percSurf;
            FitR2_percSurf(i) = baselines.R2_percSurf;
            FitN1_percSurf(i) = baselines.N1_percSurf;
            FitR1_Surf(i) = baselines.R1_surf;
            FitR2_Surf(i) = baselines.R2_surf;
            FitN1_Surf(i) = baselines.N1_surf;
            FitR1_Rab4(i) = baselines.R1_Rab4;
            FitR2_Rab4(i) = baselines.R2_Rab4;
            FitN1_Rab4(i) = baselines.N1_Rab4;
            FitR1_Rab11(i) = baselines.R1_Rab11;
            FitR2_Rab11(i) = baselines.R2_Rab11;
            FitN1_Rab11(i) = baselines.N1_Rab11;
            FitR1_Int(i) = baselines.R1_int;
            FitR2_Int(i) = baselines.R2_int;
            FitN1_Int(i) = baselines.N1_int;
            
            FitR1_2hrCHX(i) = R1_sim(4)/baselines.R1_total;
            FitR2_2hrCHX(i) = R2_sim(4)/baselines.R2_total;
            FitN1_2hrCHX(i) = N1_sim(4)/baselines.N1_total;
            
            FitR1_siRNA(i) = R1_sim(8)/baselines.R1_total;
            FitR2_siRNA(i) = R2_sim(8)/baselines.R2_total;
            FitN1_siRNA(i) = N1_sim(8)/baselines.N1_total;            
            
            disp(['simulations ' num2str((i/size(fits,2))*100) ' % complete']);
        end

        f1 = figure;
        subplot(4,3,1);
        ax1 = scatter(FitR1_percSurf,FitRunCost); set(gca,'yscale','log');
        subplot(4,3,2);
        ax2 = scatter(FitR2_percSurf,FitRunCost); set(gca,'yscale','log');
        subplot(4,3,3);
        ax3 = scatter(FitN1_percSurf,FitRunCost); set(gca,'yscale','log');
        
        subplot(4,3,4);
        ax1 = scatter(FitR1_Surf,FitRunCost); set(gca,'yscale','log');
        subplot(4,3,5);
        ax2 = scatter(FitR2_Surf,FitRunCost); set(gca,'yscale','log');
        subplot(4,3,6);
        ax3 = scatter(FitN1_Surf,FitRunCost); set(gca,'yscale','log');
        
        subplot(4,3,7);
        ax1 = scatter(FitR1_2hrCHX,FitRunCost); set(gca,'yscale','log');
        subplot(4,3,8);
        ax2 = scatter(FitR2_2hrCHX,FitRunCost); set(gca,'yscale','log');
        subplot(4,3,9);
        ax3 = scatter(FitN1_2hrCHX,FitRunCost); set(gca,'yscale','log');

        subplot(4,3,10);
        ax1 = scatter(FitR1_siRNA,FitRunCost); set(gca,'yscale','log');
        subplot(4,3,11);
        ax2 = scatter(FitR2_siRNA,FitRunCost); set(gca,'yscale','log');
        subplot(4,3,12);
        ax3 = scatter(FitN1_siRNA,FitRunCost); set(gca,'yscale','log');

        f2 = figure;
        for i = 1:15
            subplot(3,5,i);
            ax1 = scatter(squeeze(inits(i,:)),squeeze(fits(i,:)),100./sqrt(FitRunCost),'MarkerFaceColor','b','MarkerEdgeColor','b','MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2);
%             ax1 = scatter(squeeze(fits(i,:)),sqrt(FitRunCost),100./sqrt(FitRunCost),'MarkerFaceColor','b','MarkerEdgeColor','b','MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2);
            set(gca,'Xscale','log');           
            set(gca,'Yscale','log');           
            xlim([1e-6 1e-1]);
            ylim([1e-6 1e-1]);
        end

        f3 = figure;
        fitsR1 = fits([1:5 16],:);
        for i = 1:6
        for j = 1:6
            subplot(6,6,(j-1)*6+i);
            ax1 = scatter(squeeze(fitsR1(i,:)),squeeze(fitsR1(j,:)),100./sqrt(FitRunCost),'MarkerFaceColor','b','MarkerEdgeColor','b','MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2);
            set(gca,'Xscale','log');           
            set(gca,'Yscale','log');           
            xlim([1e-6 1e-1]);
            ylim([1e-6 1e-1]);
        end
        end
        R1_corr = corrcoef(fitsR1');
        
        f4 = figure;
        fitsR2 = fits([6:10 17],:);
        for i = 1:6
        for j = 1:6
            subplot(6,6,(j-1)*6+i);
            ax1 = scatter(squeeze(fitsR2(i,:)),squeeze(fitsR2(j,:)),100./sqrt(FitRunCost),'MarkerFaceColor','b','MarkerEdgeColor','b','MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2);
            set(gca,'Xscale','log');           
            set(gca,'Yscale','log');           
            xlim([1e-6 1e-1]);
            ylim([1e-6 1e-1]);
        end
        end
        R2_corr = corrcoef(fitsR2');

        f5 = figure;
        fitsN1 = fits([11:15 18],:);
        for i = 1:6
        for j = 1:6
            subplot(6,6,(j-1)*6+i);
            ax1 = scatter(squeeze(fitsN1(i,:)),squeeze(fitsN1(j,:)),100./sqrt(FitRunCost),'MarkerFaceColor','b','MarkerEdgeColor','b','MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2);
            set(gca,'Xscale','log');           
            set(gca,'Yscale','log');           
            xlim([1e-6 1e-1]);
            ylim([1e-6 1e-1]);
        end
        end
        N1_corr = corrcoef(fitsN1');
        
        paramnames = {'kR1Rab5a','kR1Rab4at7a','kR1Rab4a','kR1Rab4at11a','kR1Rab11a',...
            'kR2Rab5a','kR2Rab4at7a','kR2Rab4a','kR2Rab4at11a','kR2Rab11a',...
            'kN1Rab5a','kN1Rab4at7a','kN1Rab4a','kN1Rab4at11a','kN1Rab11a',...
            'kR1prod','kR2prod','kN1prod'}; 
        
        reframe = [1 6 11 2 7 12 3 8 13 4 9 14 5 10 15 16 17 18];
        fitsreframed = fits(reframe,:);
        paramnamesreframed = paramnames(reframe);

        fits_table = table(FitRunCost,FitR1_percSurf,FitR2_percSurf,FitN1_percSurf,...
                            FitR1_Surf,FitR2_Surf,FitN1_Surf,...
                            FitR1_2hrCHX,FitR2_2hrCHX,FitN1_2hrCHX,...
                            FitR1_siRNA,FitR2_siRNA,FitN1_siRNA);

        writetable(fits_table, "Fit_outcomes.csv");
        params_table = table(fits');
        params_table = splitvars(params_table,"Var1","NewVariableNames",paramnames);
        writetable(params_table, "Fit_params.csv");

        
        f6 = figure;
        fitsR1fit = [FitR1_percSurf FitR1_Surf FitR1_2hrCHX FitR1_siRNA];
        for i = 1:4
        for j = 1:4
            subplot(4,4,(j-1)*4+i);
            ax1 = scatter(squeeze(fitsR1fit(:,i)),squeeze(fitsR1fit(:,j)),100./sqrt(FitRunCost),'MarkerFaceColor','b','MarkerEdgeColor','b','MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2);
%             set(gca,'Xscale','log');           
%             set(gca,'Yscale','log');           
%             xlim([1e-6 1e-1]);
%             ylim([1e-6 1e-1]);
        end
        end
        R1fit_corr = corrcoef(fitsR1fit);
        
        f7 = figure;
        fitsR2fit = [FitR2_percSurf FitR2_Surf FitR2_2hrCHX FitR2_siRNA];
        for i = 1:4
        for j = 1:4
            subplot(4,4,(j-1)*4+i);
            ax1 = scatter(squeeze(fitsR2fit(:,i)),squeeze(fitsR2fit(:,j)),100./sqrt(FitRunCost),'MarkerFaceColor','b','MarkerEdgeColor','b','MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2);
%             set(gca,'Xscale','log');           
%             set(gca,'Yscale','log');           
%             xlim([1e-6 1e-1]);
%             ylim([1e-6 1e-1]);
        end
        end
        R2fit_corr = corrcoef(fitsR2fit);

        f8 = figure;
        fitsN1fit = [FitN1_percSurf FitN1_Surf FitN1_2hrCHX FitN1_siRNA];
        for i = 1:4
        for j = 1:4
            subplot(4,4,(j-1)*4+i);
            ax1 = scatter(squeeze(fitsN1fit(:,i)),squeeze(fitsN1fit(:,j)),100./sqrt(FitRunCost),'MarkerFaceColor','b','MarkerEdgeColor','b','MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2);
%             set(gca,'Xscale','log');           
%             set(gca,'Yscale','log');           
%             xlim([1e-6 1e-1]);
%             ylim([1e-6 1e-1]);
        end
        end
        N1fit_corr = corrcoef(fitsN1fit);
        
        f9 = figure;
        fitsfit = [fitsR1fit fitsR2fit fitsN1fit];
        for i = 1:12
            subplot(3,4,i);
            ax1 = scatter(squeeze(fitsfit(:,i)),FitRunCost,100./sqrt(FitRunCost),'MarkerFaceColor','b','MarkerEdgeColor','b','MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2);
%             set(gca,'Xscale','log');           
            set(gca,'Yscale','log');           
%             xlim([1e-6 1e-1]);
%             ylim([1e-6 1e-1]);
        end
        fitsfit_corr = corrcoef(fitsfit);

        levelsout = [FitR1_Surf FitR2_Surf FitN1_Surf];
        csvwrite('ManyFits_SurfaceRec.csv',levelsout);
        levelsout = [FitR1_Int FitR2_Int FitN1_Int];
        csvwrite('ManyFits_InternalRec.csv',levelsout);
        levelsout = [FitR1_Rab4 FitR2_Rab4 FitN1_Rab4];
        csvwrite('ManyFits_Rab4Rec.csv',levelsout);
        levelsout = [FitR1_Rab11 FitR2_Rab11 FitN1_Rab11];
        csvwrite('ManyFits_Rab11Rec.csv',levelsout);
        
        
    case 3
        
        % define global variables for tracking optimization loops (for debugging purposes)
        global r1loop
        global ploop
        r1loop = 0;
        ploop = 0;

% SELECT WHICH ONE PARAMETER SET TO RUN:        
%         p1 = readmatrix('LWCfit2015.csv');
%         p1 = readmatrix('SSfit2022.csv');
        p1 = readmatrix('SSfit2023.csv');

        parameters.kR1Rab5a     = p1(1); 
        parameters.kR1Rab4at7a  = p1(2); 
        parameters.kR1Rab4a     = p1(3); 
        parameters.kR1Rab4at11a = p1(4); 
        parameters.kR1Rab11a    = p1(5); 

        parameters.kR1N1Rab5a     = p1(1);
        parameters.kR1N1Rab4at7a  = p1(2); 
        parameters.kR1N1Rab4a     = p1(3); 
        parameters.kR1N1Rab4at11a = p1(4); 
        parameters.kR1N1Rab11a    = p1(5); 
        
        parameters.kR2Rab5a     = p1(6); 
        parameters.kR2Rab4at7a  = p1(7); 
        parameters.kR2Rab4a     = p1(8); 
        parameters.kR2Rab4at11a = p1(9); 
        parameters.kR2Rab11a    = p1(10);
        
        parameters.kN1Rab5a     = p1(11);
        parameters.kN1Rab4at7a  = p1(12);
        parameters.kN1Rab4a     = p1(13);
        parameters.kN1Rab4at11a = p1(14);
        parameters.kN1Rab11a    = p1(15);

        parameters.kR1prod = p1(16); 
        parameters.kR2prod = p1(17); 
        parameters.kN1prod = p1(18); 

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

% Alternative - particleswarm optimization (inefficient for this monotonic problem)
%             global PStargets
%             global PSparameters
%             global PSmodel
%             PStargets = recsurftargets;
%             PSparameters = parameters;
%             PSmodel = model;
% %             options = optimoptions('particleswarm','SwarmSize',100,'HybridFcn',@fmincon,'Display','iter');
%             options = optimoptions('particleswarm','SwarmSize',100,'Display','iter');
%             [optimalR, fval, exitflag, output] = particleswarm(@CostFxn_ReceptorProductionPS,length(p2),lb,ub,options);

            parameters.kR1prod = optimalR(1);
            parameters.kR2prod = optimalR(2);
            parameters.kN1prod = optimalR(3);
        end
        
        [timepoints, species_out, observables_out] = SimToSteadyState(parameters,model);
        speciesInit = species_out(end,:); % we just ran everything to steady state, so use the outputs as new inputs
        ConcOut=CalcOutputs(observables_out,suppressligands);
        baselines = CalcBaselines(ConcOut);

        % RUN Unliganded EXPERIMENTS and output figures
        % severity of knockdowns 
        perturbs.R1reductionByCHX = 1.00; % 0 = complete shutdown of production
        perturbs.R2reductionByCHX = 1.00; % 0.5 = half of production
        perturbs.N1reductionByCHX = 1.00; % 1.0 = no change to production

        perturbs.Rab4reductionBysiRNA    = 1.00;
        perturbs.Rab11reductionBysiRNA   = 1.00;

        testCHQ = [1 0.6 0.58 0.56 0.54 0.52];
        ratesout = [];
        VEGFR1outputs = [];
        for j = 1:length(testCHQ)
            perturbs.Rab7reductionByCHQ    = testCHQ(j);
            runparams = parameters;
            runparams.kR1Rab4at7a = parameters.kR1Rab4at7a * perturbs.Rab7reductionByCHQ;
            runparams.kR2Rab4at7a = parameters.kR2Rab4at7a * perturbs.Rab7reductionByCHQ;
            runparams.kN1Rab4at7a = parameters.kN1Rab4at7a * perturbs.Rab7reductionByCHQ;
            runparams.kR1N1Rab4at7a = parameters.kR1N1Rab4at7a * perturbs.Rab7reductionByCHQ;
        
            Timepoints = 0:60*5:3600*24;
            [~, timepoints, species_out, observables_out] = eval(append(model,"(Timepoints', speciesInit, runparams, 1)"));
            ConcOutCHQ(j)=CalcOutputs(observables_out,suppressligands);
            baselines = CalcBaselines(ConcOut);
            [R1rates(j),R2rates(j),N1rates(j)] = CalcRates(ConcOutCHQ(j),runparams,perturbs);
            ratenames = fieldnames(R1rates);
            for k = 1:length(ratenames)
                R1ratesout(k) = R1rates(j).(ratenames{k})(end);
                R2ratesout(k) = R2rates(j).(ratenames{k})(end);
                N1ratesout(k) = N1rates(j).(ratenames{k})(end);
            end
            % 4 hr timepoint = point 49; 18 hr timepoint = point 217
            VEGFR1outputs = [VEGFR1outputs [ConcOutCHQ(j).R1_surf(49)/ConcOutCHQ(j).R1_surf(1);ConcOutCHQ(j).R1_int(49)/ConcOutCHQ(j).R1_int(1);...
                            ConcOutCHQ(j).R1_surf(217)/ConcOutCHQ(j).R1_surf(1);ConcOutCHQ(j).R1_int(217)/ConcOutCHQ(j).R1_int(1)]];
            ratesout = [ratesout; [R1ratesout' R2ratesout' N1ratesout']];
        end
        csvwrite('EndpointRates.csv',ratesout);
        csvwrite('VEGFR1_4h_18h.csv',VEGFR1outputs*100);
        ConcGraphR1 = []; ConcGraphR2 = []; ConcGraphN1 = [];
        for j = 1:5; 
            ConcGraphR1 = [ConcGraphR1 ConcOutCHQ(j).R1_total'./ConcOutCHQ(1).R1_total'*100]; 
            ConcGraphR2 = [ConcGraphR2 ConcOutCHQ(j).R2_total'./ConcOutCHQ(1).R2_total'*100]; 
            ConcGraphN1 = [ConcGraphN1 ConcOutCHQ(j).N1_total'./ConcOutCHQ(1).N1_total'*100]; 
        end
        figure; ax1 = subplot(3,1,1); plot(timepoints, ConcGraphR1);
                ax2 = subplot(3,1,2); plot(timepoints, ConcGraphR2);
                ax3 = subplot(3,1,3); plot(timepoints, ConcGraphN1);
        csvwrite('ReceptorCurvesCHQ.csv',[[timepoints/3600 ConcGraphR1];[timepoints/3600 ConcGraphR2];[timepoints/3600 ConcGraphN1]]);                


    case 4
        
        % define global variables for tracking optimization loops (for debugging purposes)
        global r1loop
        global ploop
        r1loop = 0;
        ploop = 0;

% SELECT WHICH ONE PARAMETER SET TO RUN:        
%         p1 = readmatrix('LWCfit2015.csv');
%         p1 = readmatrix('SSfit2022.csv');
        p1 = readmatrix('SSfit2023.csv');

        parameters.kR1Rab5a     = p1(1); 
        parameters.kR1Rab4at7a  = p1(2); 
        parameters.kR1Rab4a     = p1(3); 
        parameters.kR1Rab4at11a = p1(4); 
        parameters.kR1Rab11a    = p1(5); 

        parameters.kR1N1Rab5a     = p1(1);
        parameters.kR1N1Rab4at7a  = p1(2); 
        parameters.kR1N1Rab4a     = p1(3); 
        parameters.kR1N1Rab4at11a = p1(4); 
        parameters.kR1N1Rab11a    = p1(5); 
        
        parameters.kR2Rab5a     = p1(6); 
        parameters.kR2Rab4at7a  = p1(7); 
        parameters.kR2Rab4a     = p1(8); 
        parameters.kR2Rab4at11a = p1(9); 
        parameters.kR2Rab11a    = p1(10);
        
        parameters.kN1Rab5a     = p1(11);
        parameters.kN1Rab4at7a  = p1(12);
        parameters.kN1Rab4a     = p1(13);
        parameters.kN1Rab4at11a = p1(14);
        parameters.kN1Rab11a    = p1(15);

        parameters.kR1prod = p1(16); 
        parameters.kR2prod = p1(17); 
        parameters.kN1prod = p1(18); 

%         coupling_modifiers = [0 logspace(-5,1,19)];
        coupling_modifiers = [0 logspace(0,4,21)];

        outputN1prod = zeros(length(coupling_modifiers),1);
        outputN1baseline = outputN1prod;
        outputN1underRabs = outputN1prod;
        outputN1underRabRatio = outputN1prod;
        outputkN1R1on = outputN1prod;
        baselinekN1R1on_surf = parameters.kN1R1on_surf;
        baselinekN1R1on_rab4 = parameters.kN1R1on_rab4;
        baselinekN1R1on_rab11 = parameters.kN1R1on_rab11;
        for i = 1:length(coupling_modifiers)
            parameters.kN1R1on_surf = baselinekN1R1on_surf*coupling_modifiers(i);
            parameters.kN1R1on_rab4 = baselinekN1R1on_rab4*coupling_modifiers(i);
            parameters.kN1R1on_rab11 = baselinekN1R1on_rab11*coupling_modifiers(i);
            outputkN1R1on(i) = parameters.kN1R1on_surf*1000;
        
%   Optimize production rates to surface receptor levels, if needed
        if optimizeR == 1
            currentR1prodguess = parameters.kR1prod;
            currentR2prodguess = parameters.kR2prod;
            currentN1prodguess = parameters.kN1prod;
            p2(1)=currentR1prodguess;
            p2(2)=currentR2prodguess;
            p2(3)=currentN1prodguess; 
            initguesses = p2;
            lb=p2*.00000001;     % lower bound
            ub=p2*10000000;     % upper bound
            options = optimoptions('lsqnonlin','Display','iter'); % display output
            optimalR = lsqnonlin(@CostFxn_ReceptorProduction,p2,lb,ub,options,recsurftargets,parameters,model);

% Alternative - particleswarm optimization (inefficient for this monotonic problem)
%             global PStargets
%             global PSparameters
%             global PSmodel
%             PStargets = recsurftargets;
%             PSparameters = parameters;
%             PSmodel = model;
% %             options = optimoptions('particleswarm','SwarmSize',100,'HybridFcn',@fmincon,'Display','iter');
%             options = optimoptions('particleswarm','SwarmSize',100,'Display','iter');
%             [optimalR, fval, exitflag, output] = particleswarm(@CostFxn_ReceptorProductionPS,length(p2),lb,ub,options);

            parameters.kR1prod = optimalR(1);
            parameters.kR2prod = optimalR(2);
            parameters.kN1prod = optimalR(3);
            outputN1prod(i) = parameters.kN1prod;
        end

    % RUN Unliganded EXPERIMENTS and output figures
        % severity of knockdowns 
        perturbs.R1reductionByCHX = 0.00; % 0 = complete shutdown of production
        perturbs.R2reductionByCHX = 0.00; % 0.5 = half of production
        perturbs.N1reductionByCHX = 0.00; % 1.0 = no change to production

        perturbs.Rab4reductionBysiRNA    = 1 - 0.80;  %80% reduction
        perturbs.Rab11reductionBysiRNA   = 1 - 0.80;

        [timepoints, species_out, observables_out] = SimToSteadyState(parameters,model);
        speciesInit = species_out(end,:); % we just ran everything to steady state, so use the outputs as new inputs
        ConcOut=CalcOutputs(observables_out,suppressligands);
        baselines = CalcBaselines(ConcOut);
        
        outputN1baseline(i) = baselines.N1_total;
        
        perturbsnone = perturbs;
        perturbfields = fieldnames(perturbs);
        for k = 1:length(perturbfields)
            perturbsnone.(perturbfields{k}) = 1;
        end
        [R1rates,R2rates,N1rates] = CalcRates(ConcOut,parameters,perturbsnone);
        ratenames = fieldnames(R1rates);
        for k = 1:length(ratenames)
            R1ratesout(k) = R1rates.(ratenames{k})(end);
            R2ratesout(k) = R2rates.(ratenames{k})(end);
            N1ratesout(k) = N1rates.(ratenames{k})(end);
        end
        ratesout = [R1ratesout' R2ratesout' N1ratesout'];
        csvwrite('EndpointRates.csv',ratesout);
        
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
        
        createFigs = 0;
        [R1_exp,R1_sim,R2_exp,R2_sim,N1_exp,N1_sim] = UnligandedHUVECExperimentsAll(parameters,speciesInit,baselines,perturbs,recsurftargets,model,createFigs);
        outputN1underRabs(i) = N1_sim(8);
        end
        
        outputN1 = [outputkN1R1on outputN1prod];
        csvwrite('NRP1_production_rates.csv',outputN1);
        
        outputN1 = [outputkN1R1on outputN1baseline];
        csvwrite('NRP1_wholecell_baselines.csv',outputN1);

        outputN1 = [outputkN1R1on outputN1underRabs];
        csvwrite('NRP1_wholecell_underRab.csv',outputN1);
        
        outputN1underRabRatio = outputN1underRabs./outputN1baseline;
        outputN1 = [outputkN1R1on outputN1underRabRatio];
        csvwrite('NRP1_wholecell_underRab_ratio.csv',outputN1);        

        
    case 5
        
        % define global variables for tracking optimization loops (for debugging purposes)
        global r1loop
        global ploop
        r1loop = 0;
        ploop = 0;

% SELECT WHICH ONE PARAMETER SET TO RUN:        
%         p1 = readmatrix('LWCfit2015.csv');
%         p1 = readmatrix('SSfit2022.csv');
        p1 = readmatrix('SSfit2023.csv');

        parameters.kR1Rab5a     = p1(1); 
        parameters.kR1Rab4at7a  = p1(2); 
        parameters.kR1Rab4a     = p1(3); 
        parameters.kR1Rab4at11a = p1(4); 
        parameters.kR1Rab11a    = p1(5); 

        parameters.kR1N1Rab5a     = p1(1);
        parameters.kR1N1Rab4at7a  = p1(2); 
        parameters.kR1N1Rab4a     = p1(3); 
        parameters.kR1N1Rab4at11a = p1(4); 
        parameters.kR1N1Rab11a    = p1(5); 
        
        parameters.kR2Rab5a     = p1(6); 
        parameters.kR2Rab4at7a  = p1(7); 
        parameters.kR2Rab4a     = p1(8); 
        parameters.kR2Rab4at11a = p1(9); 
        parameters.kR2Rab11a    = p1(10);
        
        parameters.kN1Rab5a     = p1(11);
        parameters.kN1Rab4at7a  = p1(12);
        parameters.kN1Rab4a     = p1(13);
        parameters.kN1Rab4at11a = p1(14);
        parameters.kN1Rab11a    = p1(15);

        parameters.kR1prod = p1(16); 
        parameters.kR2prod = p1(17); 
        parameters.kN1prod = p1(18); 

%         varyparam = "R1R1";
%         varyparam = "R2R2";
        varyparam = "N1R1";
        
%         parameters.kN1R1on  = 0; 
        
        if optimizeR == 0 % modifiers below are used to check on stoichiometry and balance of the model
            % only implement when *not* re-optimizing receptor production rates

            % Production & degradation - set to zero to check static equilibrium
            
            % no degradation if following are set to zero
            parameters.kR1Rab4at7a = 0;
            parameters.kR2Rab4at7a = 0;
            parameters.kN1Rab4at7a = 0;
            parameters.kR1N1Rab4at7a = 0; 

            % no new production if following are set to zero
            parameters.kR1prod = 0; 
            parameters.kR2prod = 0; 
            parameters.kN1prod = 0; 

%             % dimerization - set to zero to switch off R1 dimerization 
%             parameters.kR1R1on = 0; 

            % Trafficking - set to ero to switch off trafficking processes
            
%             % internalization - set to zero to keep everything on the surface 
%             parameters.kR1Rab5a = 0;
%             parameters.kR2Rab5a = 0;
%             parameters.kN1Rab5a = 0;
%             parameters.kR1N1Rab5a = parameters.kR1Rab5a;
            
%             % Rab4 recycling - set to zero to switch off process
%             parameters.kR1Rab4a = 0;
%             parameters.kR2Rab4a = 0;
%             parameters.kN1Rab4a = 0;
%             parameters.kR1N1Rab4a = parameters.kR1Rab4a; 
% 
%             % Rab4 to Rab11 transfer - set to zero to switch off process
%             parameters.kR1Rab4at11a = 0;
%             parameters.kR2Rab4at11a = 0;
%             parameters.kN1Rab4at11a = 0;
%             parameters.kR1N1Rab4at11a = parameters.kR1Rab4at11a; 
%  
%             % Rab11 recycling - set to zero to switch off process
%             parameters.kR1Rab11a = 0;
%             parameters.kR2Rab11a = 0;
%             parameters.kN1Rab11a = 0;
%             parameters.kR1N1Rab11a = parameters.kR1Rab11a; 

        end

%         coupling_modifiers = [logspace(2,-2,21) 0];
%         coupling_modifiers = [logspace(1,-1,11) 0];
        coupling_modifiers = [0 logspace(-2,2,21)];

        outputR1prod = zeros(length(coupling_modifiers),1);
        outputR1tot_surf =  outputR1prod;
        outputR1tot_Rab4 =  outputR1prod;
        outputR1tot_Rab11 =  outputR1prod;
        outputR1dimers = outputR1prod;
        outputR1dimeric = outputR1prod;
        outputR1dimeric_surf = outputR1prod;
        outputR1dimeric_Rab4 = outputR1prod;
        outputR1dimeric_Rab11 = outputR1prod;
        outputkR1R1on = outputR1prod;
        baselinekR1R1on_surf = parameters.kR1R1on_surf;
        baselinekR1R1on_rab4 = parameters.kR1R1on_rab4;
        baselinekR1R1on_rab11 = parameters.kR1R1on_rab11;

        outputR2prod = outputR1prod;
        outputR2tot_surf =  outputR1prod;
        outputR2tot_Rab4 =  outputR1prod;
        outputR2tot_Rab11 =  outputR1prod;
        outputR2dimers = outputR1prod;
        outputR2dimeric = outputR1prod;
        outputR2dimeric_surf = outputR1prod;
        outputR2dimeric_Rab4 = outputR1prod;
        outputR2dimeric_Rab11 = outputR1prod;
        outputkR2R2on = outputR1prod;
        baselinekR2R2on_surf = parameters.kR2R2on_surf;
        baselinekR2R2on_rab4 = parameters.kR2R2on_rab4;
        baselinekR2R2on_rab11 = parameters.kR2R2on_rab11;

        outputN1prod = outputR1prod;
        outputN1tot_surf =  outputR1prod;
        outputN1tot_Rab4 =  outputR1prod;
        outputN1tot_Rab11 =  outputR1prod;
        outputR1_R1N1fraction = outputR1prod;
        outputN1_R1N1fraction = outputR1prod;
        outputR1_R1N1fraction_surf = outputR1prod;
        outputR1_R1N1fraction_Rab4 = outputR1prod;
        outputR1_R1N1fraction_Rab11 = outputR1prod;
        outputN1_R1N1fraction_surf = outputR1prod;
        outputN1_R1N1fraction_Rab4 = outputR1prod;
        outputN1_R1N1fraction_Rab11 = outputR1prod;

        outputR1_R1N1total = outputR1prod;
        outputN1_R1N1total = outputR1prod;
        outputR1_R1N1surf = outputR1prod;
        outputR1_R1N1Rab4 = outputR1prod;
        outputR1_R1N1Rab11 = outputR1prod;
        outputN1_R1N1surf = outputR1prod;
        outputN1_R1N1Rab4 = outputR1prod;
        outputN1_R1N1Rab11 = outputR1prod;

        outputkN1R1on = outputR1prod;
        baselinekN1R1on_surf = parameters.kN1R1on_surf;
        baselinekN1R1on_rab4 = parameters.kN1R1on_rab4;
        baselinekN1R1on_rab11 = parameters.kN1R1on_rab11;
        
        for i = 1:length(coupling_modifiers)

            switch varyparam
                case "R1R1"
                    parameters.kR1R1on_surf = baselinekR1R1on_surf*coupling_modifiers(i);
                    parameters.kR1R1on_rab4 = baselinekR1R1on_rab4*coupling_modifiers(i);
                    parameters.kR1R1on_rab11 = baselinekR1R1on_rab11*coupling_modifiers(i);
                case "R2R2"
                    parameters.kR2R2on_surf = baselinekR2R2on_surf*coupling_modifiers(i);
                    parameters.kR2R2on_rab4 = baselinekR2R2on_rab4*coupling_modifiers(i);
                    parameters.kR2R2on_rab11 = baselinekR2R2on_rab11*coupling_modifiers(i);
                case "N1R1"
                    parameters.kN1R1on_surf = baselinekN1R1on_surf*coupling_modifiers(i);
                    parameters.kN1R1on_rab4 = baselinekN1R1on_rab4*coupling_modifiers(i);
                    parameters.kN1R1on_rab11 = baselinekN1R1on_rab11*coupling_modifiers(i);
            end

            outputkR1R1on(i) = parameters.kR1R1on_surf*1000;
            outputkR2R2on(i) = parameters.kR2R2on_surf*1000;
            outputkN1R1on(i) = parameters.kN1R1on_surf*1000;
        
%   Optimize production rates to surface receptor levels, if needed
        if optimizeR == 1
            currentR1prodguess = parameters.kR1prod;
            currentR2prodguess = parameters.kR2prod;
            currentN1prodguess = parameters.kN1prod;
            p2(1)=currentR1prodguess;
            p2(2)=currentR2prodguess;
            p2(3)=currentN1prodguess; 
            initguesses = p2;
            lb=p2*.00000001;     % lower bound
            ub=p2*10000000;     % upper bound
            options = optimoptions('lsqnonlin','Display','iter'); % display output
            optimalR = lsqnonlin(@CostFxn_ReceptorProduction,p2,lb,ub,options,recsurftargets,parameters,model);

            parameters.kR1prod = optimalR(1);
            parameters.kR2prod = optimalR(2);
            parameters.kN1prod = optimalR(3);
            outputR1prod(i) = parameters.kR1prod;
            outputR2prod(i) = parameters.kR2prod;
            outputN1prod(i) = parameters.kN1prod;
        end

    % RUN Unliganded EXPERIMENTS and output figures
        % severity of knockdowns 
        perturbs.R1reductionByCHX = 0.00; % 0 = complete shutdown of production
        perturbs.R2reductionByCHX = 0.00; % 0.5 = half of production
        perturbs.N1reductionByCHX = 0.00; % 1.0 = no change to production

        perturbs.Rab4reductionBysiRNA    = 1 - 0.80;  %80% reduction
        perturbs.Rab11reductionBysiRNA   = 1 - 0.80;

        [timepoints, species_out, observables_out] = SimToSteadyState(parameters,model);
        speciesInit = species_out(end,:); % we just ran everything to steady state, so use the outputs as new inputs
        ConcOut=CalcOutputs(observables_out,suppressligands);
        baselines = CalcBaselines(ConcOut);
        
        outputR1tot_surf(i) = baselines.R1_surf;
        outputR1tot_Rab4(i) = baselines.R1_Rab4;
        outputR1tot_Rab11(i) = baselines.R1_Rab11;
        outputR2tot_surf(i) = baselines.R2_surf;
        outputR2tot_Rab4(i) = baselines.R2_Rab4;
        outputR2tot_Rab11(i) = baselines.R2_Rab11;
        outputN1tot_surf(i) = baselines.N1_surf;
        outputN1tot_Rab4(i) = baselines.N1_Rab4;
        outputN1tot_Rab11(i) = baselines.N1_Rab11;

        outputR1dimeric(i) = baselines.R1_dimericfraction;
        outputR1dimeric_surf(i) = baselines.R1_dimericfraction_surf;
        outputR1dimeric_Rab4(i) = baselines.R1_dimericfraction_Rab4;
        outputR1dimeric_Rab11(i) = baselines.R1_dimericfraction_Rab11;
        outputR2dimeric(i) = baselines.R2_dimericfraction;
        outputR2dimeric_surf(i) = baselines.R2_dimericfraction_surf;
        outputR2dimeric_Rab4(i) = baselines.R2_dimericfraction_Rab4;
        outputR2dimeric_Rab11(i) = baselines.R2_dimericfraction_Rab11;
        outputR1dimers(i) = baselines.R1_dimers;
        outputR2dimers(i) = baselines.R2_dimers;
        outputR1_R1N1fraction(i) = baselines.R1_fractioninR1N1;
        outputR1_R1N1fraction_surf(i) = baselines.R1_fractioninR1N1_surf;
        outputR1_R1N1fraction_Rab4(i) = baselines.R1_fractioninR1N1_Rab4;
        outputR1_R1N1fraction_Rab11(i) = baselines.R1_fractioninR1N1_Rab11;
        outputN1_R1N1fraction(i) = baselines.N1_fractioninR1N1;
        outputN1_R1N1fraction_surf(i) = baselines.N1_fractioninR1N1_surf;
        outputN1_R1N1fraction_Rab4(i) = baselines.N1_fractioninR1N1_Rab4;
        outputN1_R1N1fraction_Rab11(i) = baselines.N1_fractioninR1N1_Rab11;

        outputR1_R1N1total(i) = baselines.R1_in_R1N1_total;
        outputR1_R1N1surf(i) = baselines.R1_in_R1N1_surf;
        outputR1_R1N1Rab4(i) = baselines.R1_in_R1N1_Rab4;
        outputR1_R1N1Rab11(i) = baselines.R1_in_R1N1_Rab11;
        outputN1_R1N1total(i) = baselines.N1_in_R1N1_total;
        outputN1_R1N1surf(i) = baselines.N1_in_R1N1_surf;
        outputN1_R1N1Rab4(i) = baselines.N1_in_R1N1_Rab4;
        outputN1_R1N1Rab11(i) = baselines.N1_in_R1N1_Rab11;

        perturbsnone = perturbs;
        perturbfields = fieldnames(perturbs);
        for k = 1:length(perturbfields)
            perturbsnone.(perturbfields{k}) = 1;
        end
        [R1rates,R2rates,N1rates] = CalcRates(ConcOut,parameters,perturbsnone);
        ratenames = fieldnames(R1rates);
        for k = 1:length(ratenames)
            R1ratesout(k) = R1rates.(ratenames{k})(end);
            R2ratesout(k) = R2rates.(ratenames{k})(end);
            N1ratesout(k) = N1rates.(ratenames{k})(end);
        end
        ratesout = [R1ratesout' R2ratesout' N1ratesout'];
        csvwrite('EndpointRates.csv',ratesout);
        
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
                
        end

        outputR1tot = outputR1tot_surf+outputR1tot_Rab4+outputR1tot_Rab11;
        outputR2tot = outputR2tot_surf+outputR2tot_Rab4+outputR2tot_Rab11;
        outputN1tot = outputN1tot_surf+outputN1tot_Rab4+outputN1tot_Rab11;

        switch varyparam
            case "R1R1"
                varyingrate = outputkR1R1on;
            case "R2R2"
                varyingrate = outputkR2R2on;
            case "N1R1"
                varyingrate = outputkN1R1on;
        end

        outputR1 = [varyingrate outputR1prod];
        csvwrite('Sim5_VEGFR1_production_rates.csv',outputR1);
        outputR2 = [varyingrate outputR2prod];
        csvwrite('Sim5_VEGFR2_production_rates.csv',outputR2);
        outputN1 = [varyingrate outputN1prod];
        csvwrite('Sim5_NRP1_production_rates.csv',outputN1);
        
%         outputR1 = [varyingrate outputR1dimeric];
%         csvwrite('Sim5_VEGFR1_dimericfrac_wholecell.csv',outputR1);
%         outputR1 = [varyingrate outputR1dimeric_surf];
%         csvwrite('Sim5_VEGFR1_dimericfrac_surf.csv',outputR1);
%         outputR1 = [varyingrate outputR1dimeric_Rab4];
%         csvwrite('Sim5_VEGFR1_dimericfrac_Rab4.csv',outputR1);
%         outputR1 = [varyingrate outputR1dimeric_Rab11];
%         csvwrite('Sim5_VEGFR1_dimericfrac_Rab11.csv',outputR1);
        outputR1 = [varyingrate outputR1dimeric outputR1dimeric_surf outputR1dimeric_Rab4 outputR1dimeric_Rab11];
        csvwrite('Sim5_VEGFR1_dimericfrac_bylocation.csv',outputR1);

%         outputR2 = [varyingrate outputR2dimeric];
%         csvwrite('Sim5_VEGFR2_dimericfrac_wholecell.csv',outputR2);
%         outputR2 = [varyingrate outputR2dimeric_surf];
%         csvwrite('Sim5_VEGFR2_dimericfrac_surf.csv',outputR2);
%         outputR2 = [varyingrate outputR2dimeric_Rab4];
%         csvwrite('Sim5_VEGFR2_dimericfrac_Rab4.csv',outputR2);
%         outputR2 = [varyingrate outputR2dimeric_Rab11];
%         csvwrite('Sim5_VEGFR2_dimericfrac_Rab11.csv',outputR2);
        outputR2 = [varyingrate outputR2dimeric outputR2dimeric_surf outputR2dimeric_Rab4 outputR2dimeric_Rab11];
        csvwrite('Sim5_VEGFR2_dimericfrac_bylocation.csv',outputR2);

        outputR1 = [varyingrate outputR1tot outputR1tot_surf outputR1tot_Rab4 outputR1tot_Rab11 outputR1tot_surf./outputR1tot];
        csvwrite('Sim5_VEGFR1_bylocation.csv',outputR1);
        outputR2 = [varyingrate outputR2tot outputR2tot_surf outputR2tot_Rab4 outputR2tot_Rab11 outputR2tot_surf./outputR2tot];
        csvwrite('Sim5_VEGFR2_bylocation.csv',outputR2);
        outputN1 = [varyingrate outputN1tot outputN1tot_surf outputN1tot_Rab4 outputN1tot_Rab11 outputN1tot_surf./outputN1tot];
        csvwrite('Sim5_NRP1_bylocation.csv',outputN1);

        outputR1 = [varyingrate outputR1dimers];
        csvwrite('Sim5_VEGFR1_dimers.csv',outputR1);
        outputR2 = [varyingrate outputR2dimers];
        csvwrite('Sim5_VEGFR2_dimers.csv',outputR2);

%         outputR1 = [varyingrate outputR1_R1N1fraction];
%         csvwrite('Sim5_VEGFR1_in_R1N1_frac_wholecell.csv',outputR1);
%         outputR1 = [varyingrate outputR1_R1N1fraction_surf];
%         csvwrite('Sim5_VEGFR1_in_R1N1_frac_surf.csv',outputR1);
%         outputR1 = [varyingrate outputR1_R1N1fraction_Rab4];
%         csvwrite('Sim5_VEGFR1_in_R1N1_frac_Rab4.csv',outputR1);
%         outputR1 = [varyingrate outputR1_R1N1fraction_Rab11];
%         csvwrite('Sim5_VEGFR1_in_R1N1_frac_Rab11.csv',outputR1);
        outputR1 = [varyingrate outputR1_R1N1fraction outputR1_R1N1fraction_surf outputR1_R1N1fraction_Rab4 outputR1_R1N1fraction_Rab11];
        csvwrite('Sim5_VEGFR1_in_R1N1_frac_bylocation.csv',outputR1);

%         outputN1 = [varyingrate outputN1_R1N1fraction];
%         csvwrite('Sim5_NRP1_in_R1N1_frac_wholecell.csv',outputN1);
%         outputN1 = [varyingrate outputN1_R1N1fraction_surf];
%         csvwrite('Sim5_NRP1_in_R1N1_frac_surf.csv',outputN1);
%         outputN1 = [varyingrate outputN1_R1N1fraction_Rab4];
%         csvwrite('Sim5_NRP1_in_R1N1_frac_Rab4.csv',outputN1);
%         outputN1 = [varyingrate outputN1_R1N1fraction_Rab11];
%         csvwrite('Sim5_NRP1_in_R1N1_frac_Rab11.csv',outputN1);
        outputN1 = [varyingrate outputN1_R1N1fraction outputN1_R1N1fraction_surf outputN1_R1N1fraction_Rab4 outputN1_R1N1fraction_Rab11];
        csvwrite('Sim5_NRP1_in_R1N1_frac_bylocation.csv',outputN1);

%         outputR1 = [varyingrate outputR1_R1N1total ];
%         csvwrite('Sim5_VEGFR1_in_R1N1_wholecell.csv',outputR1);
%         outputR1 = [varyingrate outputR1_R1N1surf ];
%         csvwrite('Sim5_VEGFR1_in_R1N1_surf.csv',outputR1);
%         outputR1 = [varyingrate outputR1_R1N1Rab4 ];
%         csvwrite('Sim5_VEGFR1_in_R1N1_Rab4.csv',outputR1);
%         outputR1 = [varyingrate outputR1_R1N1Rab11 ];
%         csvwrite('Sim5_VEGFR1_in_R1N1_Rab11.csv',outputR1);
        outputR1 = [varyingrate outputR1_R1N1total outputR1_R1N1surf outputR1_R1N1Rab4 outputR1_R1N1Rab11];
        csvwrite('Sim5_VEGFR1_in_R1N1_bylocation.csv',outputR1);
        
%         outputN1 = [varyingrate outputN1_R1N1total ];
%         csvwrite('Sim5_NRP1_in_R1N1_wholecell.csv',outputN1);
%         outputN1 = [varyingrate outputN1_R1N1surf ];
%         csvwrite('Sim5_NRP1_in_R1N1_surf.csv',outputN1);
%         outputN1 = [varyingrate outputN1_R1N1Rab4 ];
%         csvwrite('Sim5_NRP1_in_R1N1_Rab4.csv',outputN1);
%         outputN1 = [varyingrate outputN1_R1N1Rab11 ];
%         csvwrite('Sim5_NRP1_in_R1N1_Rab11.csv',outputN1);
        outputN1 = [varyingrate outputN1_R1N1total outputN1_R1N1surf outputN1_R1N1Rab4 outputN1_R1N1Rab11];
        csvwrite('Sim5_NRP1_in_R1N1_bylocation.csv',outputN1);
         
end

