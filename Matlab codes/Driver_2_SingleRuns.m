clear all;
close all;

% Driver for single-runs for given simulation parameters 
% VEGFR1, VEGFR2, NRP1 Trafficking Model 2022-2023
% Rab7a is the degradation compartment

%% KEY INITIALIZATIONS and OPTIONS

model = "Unligated_VEGFR_model_20230321"; % which compiled model to use
baseparams = "SSfit2022.csv"; % which baseline unliganded receptor parameters to use
parameters = base_parameters(model,baseparams); % Initialize model parameters
% parameters.kN1R1on = 0; % if preventing R1-N1 coupling (e.g. Fig S2)
timestp_stst = 3600*3; % steady-state - only output every 3 hours
% time step size for the ligand simulations is defined in the ligand simulation code (RunOneVEGFPlGFsim.m) 

% Target values for surface receptor levels
R1surf_target = 1.8;
R2surf_target = 4.9;
N1surf_target = 68;
recsurftargets = [R1surf_target; R2surf_target; N1surf_target];

suppressligands = 1; % for unliganded simulations, set = 1 to suppress ligand outputs
% for liganded simulations, set = 0 to get both

RunSim = 1;
optimizeR = 1;
% Simulations: 
% 1 = single run to steady state PLUS unliganded trafficking experiments for a given set of parameters
% 2 = multiple runs to steady state PLUS unliganded trafficking experiments
% for a given set of parameters - to calculate costs
% 3 = single run to steady state PLUS CHQ simulations
% 4 = set of runs for different R1-N1 coupling rates

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
        p1 = readmatrix('SSfit2022.csv');

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

        fits = readmatrix('75fits_20220831.csv');
        inits = readmatrix('75fits_20220831_init.csv');
        FitRunCost = zeros(size(fits,2),1);
        FitR1_percSurf = zeros(size(fits,2),1);
        FitR2_percSurf = zeros(size(fits,2),1);
        FitN1_percSurf = zeros(size(fits,2),1);

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

            % Emphasize the biotinylayion surface-internal simulations
            j=9;
            CostOut(j) = CostOut(j)*10; % give these equal weight to the other readouts; basically, try to enforce low surface R1, 17 is the number of data points in Bouchert et al (2017) paper
            CostOut(j+length(R1_sim)) = CostOut(j+length(R1_sim))*10; % give these equal weight to the other readouts; basically, try to enforce low surface R1
            CostOut(j+length(R1_sim)+length(R2_sim)) = CostOut(j+length(R1_sim)+length(R2_sim))*10; % give these equal weight to the other readouts; basically, try to enforce low surface R1

            % Emphasize the siRNA simulations
            j=6:8;
            CostOut(j) = CostOut(j)*3; 
            CostOut(j+length(R1_sim)) = CostOut(j+length(R1_sim))*3; 
            CostOut(j+length(R1_sim)+length(R2_sim)) = CostOut(j+length(R1_sim)+length(R2_sim))*3; 

            
%             FitRunCost(i) = sqrt(sum(CostOut.*CostOut));
            FitRunCost(i) = (sum(CostOut.*CostOut));
            FitR1_percSurf(i) = baselines.R1_percSurf;
            FitR2_percSurf(i) = baselines.R2_percSurf;
            FitN1_percSurf(i) = baselines.N1_percSurf;

            FitR1_2hrCHX(i) = R1_sim(4)/baselines.R1_total;
            FitR2_2hrCHX(i) = R2_sim(4)/baselines.R2_total;
            FitN1_2hrCHX(i) = N1_sim(4)/baselines.N1_total;
            
            FitR1_siRNA(i) = R1_sim(8)/baselines.R1_total;
            FitR2_siRNA(i) = R2_sim(8)/baselines.R2_total;
            FitN1_siRNA(i) = N1_sim(8)/baselines.N1_total;            
            
            disp(['simulations ' num2str((i/size(fits,2))*100) ' % complete']);
        end

        f1 = figure;
        subplot(3,3,1);
        ax1 = scatter(FitR1_percSurf,FitRunCost);
        subplot(3,3,2);
        ax2 = scatter(FitR2_percSurf,FitRunCost);
        subplot(3,3,3);
        ax3 = scatter(FitN1_percSurf,FitRunCost);

        subplot(3,3,4);
        ax1 = scatter(FitR1_2hrCHX,FitRunCost);
        subplot(3,3,5);
        ax2 = scatter(FitR2_2hrCHX,FitRunCost);
        subplot(3,3,6);
        ax3 = scatter(FitN1_2hrCHX,FitRunCost);

        subplot(3,3,7);
        ax1 = scatter(FitR1_siRNA,FitRunCost);
        subplot(3,3,8);
        ax2 = scatter(FitR2_siRNA,FitRunCost);
        subplot(3,3,9);
        ax3 = scatter(FitN1_siRNA,FitRunCost);

        f2 = figure;
        for i = 1:15
            subplot(3,5,i);
            ax1 = scatter(squeeze(inits(i,:)),squeeze(fits(i,:)),100./sqrt(FitRunCost),'MarkerFaceColor','b','MarkerEdgeColor','b','MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2);
%           ax1 = scatter(squeeze(fits(i,:)),sqrt(FitRunCost),100./sqrt(FitRunCost),'MarkerFaceColor','b','MarkerEdgeColor','b','MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2);
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

        
    case 3
        
        % define global variables for tracking optimization loops (for debugging purposes)
        global r1loop
        global ploop
        r1loop = 0;
        ploop = 0;

% SELECT WHICH ONE PARAMETER SET TO RUN:        
%         p1 = readmatrix('LWCfit2015.csv');
        p1 = readmatrix('SSfit2022.csv');

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
        for j = 1:5
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
        p1 = readmatrix('SSfit2022.csv');

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

        coupling_modifiers = [logspace(1,-5,19) 0];
        outputN1prod = zeros(length(coupling_modifiers),1);
        outputN1baseline = outputN1prod;
        outputN1underRabs = outputN1prod;
        outputN1underRabRatio = outputN1prod;
        outputkN1R1on = outputN1prod;
        baselinekN1R1on = parameters.kN1R1on;
        for i = 1:length(coupling_modifiers)
            parameters.kN1R1on = baselinekN1R1on*coupling_modifiers(i);
            outputkN1R1on(i) = parameters.kN1R1on;
        
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
end

