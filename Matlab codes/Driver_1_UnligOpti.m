clear all;
close all;

% Driver for optimization of unliganded trafficking parameters 
%  (kint, kdeg, krec4, k4to11, krec11)
%  for each of VEGFR1, VEGFR2, NRP1
%  Rab7a is the degradation compartment

%% KEY INITIALIZATIONS and OPTIONS

model = "Unligated_VEGFR_model_20230321"; % which compiled model to use
baseparams = "SSfit2022.csv"; % which baseline unliganded receptor parameters to use
parameters = base_parameters(model,baseparams); % Initialize model parameters
timestp_stst = 3600*3; % steady-state - only output every 3 hours
% time step size for the ligand simulations is defined in the ligand simulation code (RunOneVEGFPlGFsim.m) 

% Target values for surface receptor levels
R1surf_target = 1.8;
R2surf_target = 4.9;
N1surf_target = 68;
recsurftargets = [R1surf_target; R2surf_target; N1surf_target];

suppressligands = 1; % for unliganded simulations, set = 1 to suppress ligand outputs
% for liganded simulations, set = 0 to get both

% Optimization settings
InitialGuessType = 1; % 1 = starting in 2-order neighborhood of prev optimal; 0 = full five-order range
Number_Of_Optimizations = 10; % how many runs of the optimization loop to do (how many initial guesses)
rng(1); % Seed for random number generator
NumberToBurn = 0; % Number of RNG sets to burn for this seed (because previously simulated) 
Number_Of_OptimizingParams = 15;

% Initialize matrix sizes
InitGuesses = zeros(Number_Of_OptimizingParams,Number_Of_Optimizations);
FinalOptimizations = zeros(Number_Of_OptimizingParams+3,Number_Of_Optimizations); % add three for the receptor production rates

%% BURN GUESSES WHEN RE-USING SEED

if NumberToBurn > 0
    for OptiLoop = 1:NumberToBurn         
        p1 = 10.^((rand(Number_Of_OptimizingParams,1)*5)-6); % initial guesses for all optimizing parameters from 1e-6 to 1e-1
        r1 = 10.^((rand(Number_Of_OptimizingParams,1)*2)-1); % initial guesses for all optimizing parameters with values in x10 - /10 neighborhood
    end
end

%% INITIAL GUESSES

for OptiLoop = 1:Number_Of_Optimizations         
    
    p1 = 10.^((rand(Number_Of_OptimizingParams,1)*5)-6); % initial guesses for all optimizing parameters from 1e-6 to 1e-1
    r1 = 10.^((rand(Number_Of_OptimizingParams,1)*2)-1); % use this if guesses are in a narrower range

    if InitialGuessType == 1
%       Instead of full random range, here make initial guesses for all
%       optimizing parameters within one order of magnitude of guess above

%       first, load guess values
        p0 = readmatrix('SSfit2022.csv'); % 2022 parameter values (based on HUVECs)
%         p1 = readmatrix('LWCfit2015.csv'); % 2015 parameter values (based on PAECs)
        p1 = p0(1:15).*r1; % randomize initial guesses
        p1(p1>1e-1) = 1e-1; % upper limit
        p1(p1<1e-6) = 1e-6; % lower limit   
    end
    
    InitGuesses(:,OptiLoop) = p1; 

    lb=(p1*0)+(1E-6); % lower bound for optimization
    ub=(p1*0)+(1E-1); % upper bound for optimization
    
%% GLOBAL VARIABLES
    
    % define global variables for tracking optimization loops (for debugging purposes)
    global r1loop
    global ploop
    r1loop = 0;
    ploop = 0;

    % define global variables for receptor production (to keep track of recent best fits)
    global currentR1prodguess
    global currentR2prodguess
    global currentN1prodguess
    currentR1prodguess = 0.004123; %1.73e-4; %2e-3;
    currentR2prodguess = 0.001058; %1.79e-5; %2e-3;
    currentN1prodguess = 0.004121; %1.71e-3; %3.3e-3;
    if InitialGuessType == 1
        currentR1prodguess = p0(16);
        currentR2prodguess = p0(17);
        currentN1prodguess = p0(18);        
    end
    parameters.kR1prod = currentR1prodguess;
    parameters.kR2prod = currentR2prodguess;
    parameters.kN1prod = currentN1prodguess;
   
    % define a global figure - refreshed each loop to show progress to fitting parameters
    global h1
    global CheckStStateAfterRProdOpti
    global CurrGuessFig
    global ig
    CheckStStateAfterRProdOpti = figure;
    h1 = figure;
    CurrGuessFig = figure;
    figure(CurrGuessFig);
        scatter(1:15,p1,'k.');
        hold on;
        scatter(1:15,p1,'bd');
        ylim([1e-6 1e-1]);
        set(gca,'YScale','log');
%         drawnow;
    ig = p1;

%% severity of knockdowns (can put these in function call to see the effect of varying effect size)
    perturbs.R1reductionByCHX = 0.00; % 0 = complete shutdown of production
    perturbs.R2reductionByCHX = 0.00; % 0.5 = half of production
    perturbs.N1reductionByCHX = 0.00; % 1.0 = no change to production

    perturbs.Rab4reductionBysiRNA    = 1 - 0.80;  %80% reduction
    perturbs.Rab11reductionBysiRNA   = 1 - 0.80;

%% Run optimization of unliganded parameters
    options = optimoptions('lsqnonlin','Display','iter'); % display output
    optimal = lsqnonlin(@CostFxn_unligandedParams,p1,lb,ub,options, parameters, model, recsurftargets, perturbs); 
    % optimal = optimal values of p1
    saveas(h1,"Fig0_OverviewOfFits.png")

%% Assign new optimized values of parameters to the parameter set
    parameters.kR1Rab5a     = optimal(1);
    parameters.kR1Rab4at7a  = optimal(2);
    parameters.kR1Rab4a     = optimal(3);
    parameters.kR1Rab4at11a = optimal(4);
    parameters.kR1Rab11a    = optimal(5);

    parameters.kR1N1Rab5a     = optimal(1);
    parameters.kR1N1Rab4at7a  = optimal(2); 
    parameters.kR1N1Rab4a     = optimal(3); 
    parameters.kR1N1Rab4at11a = optimal(4); 
    parameters.kR1N1Rab11a    = optimal(5); 

    parameters.kR2Rab5a     = optimal(6);         
    parameters.kR2Rab4at7a  = optimal(7);     
    parameters.kR2Rab4a     = optimal(8);    
    parameters.kR2Rab4at11a = optimal(9);     
    parameters.kR2Rab11a    = optimal(10);

    parameters.kN1Rab5a     = optimal(11);         
    parameters.kN1Rab4at7a  = optimal(12);     
    parameters.kN1Rab4a     = optimal(13);    
    parameters.kN1Rab4at11a = optimal(14);     
    parameters.kN1Rab11a    = optimal(15); 

%% FINAL RUN WITH OPTIMIZED PARAMETERS

    % after optimizing the trafficking parameters, 
    % optimize the receptor production rates one last time
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

    pout = CalcPout(parameters,model,timestp_stst,suppressligands); % basically a list of key values

    % RUN KQ Unliganded EXPERIMENTS and output figures
    [timepoints, species_out, observables_out] = SimToSteadyState(parameters,model);
    speciesInit = species_out(end,:); % we just ran everything to steady state, so use the outputs as new inputs
    ConcOut=CalcOutputs(observables_out,suppressligands);
    baselines = CalcBaselines(ConcOut);
    createFigs = 1;
    [R1_exp,R1_sim,R2_exp,R2_sim,N1_exp,N1_sim] = UnligandedHUVECExperimentsAll(parameters,speciesInit,baselines,perturbs,recsurftargets,model,createFigs);

    optimalALL = [optimal;optimalR'];
    FinalOptimizations(:,OptiLoop) = optimalALL;

    save('InitialGuesses_and_OptimizedValues.mat','InitGuesses','FinalOptimizations');
    csvwrite('InitialGuessValues.csv',InitGuesses);
    csvwrite('OptimizedValues.csv',FinalOptimizations);

end


