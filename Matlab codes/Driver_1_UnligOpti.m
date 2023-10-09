clear all;
close all;

% Driver for optimization of unliganded trafficking parameters 
%  (kint, kdeg, krec4, k4to11, krec11)
%  for each of VEGFR1, VEGFR2, NRP1

if isempty(getenv('SLURM_JOB_ID'))
    poolobj = parpool(2);
else
    % If on Rockfish, use the number of CPUs requested
    poolobj = parpool(str2num(getenv('SLURM_CPUS_PER_TASK')));
end

%% KEY INITIALIZATIONS and OPTIONS

model = "Unligated_VEGFR_model_20230831"; % which compiled model to use
baseparams = "SSfit2022.csv"; % which baseline unliganded receptor parameters to use
parameters = base_parameters(model,baseparams); % Initialize model parameters
timestp_stst = 3600*3; % steady-state - only output every 3 hours

% Target values for surface receptor levels (#/cell)
R1surf_target = 1800;
R2surf_target = 4900;
N1surf_target = 68000;
recsurftargets = [R1surf_target; R2surf_target; N1surf_target];

% Severity of knockdowns
perturbs.R1reductionByCHX = 0.00; % 0 = complete shutdown of production
perturbs.R2reductionByCHX = 0.00; % 0.5 = half of production
perturbs.N1reductionByCHX = 0.00; % 1.0 = no change to production

perturbs.Rab4reductionBysiRNA    = 1 - 0.80;  %80% reduction
perturbs.Rab11reductionBysiRNA   = 1 - 0.80;
    
% Optimization settings
InitialGuessType = 1; % 1 = starting in 2-order neighborhood of prev optimal; 0 = full five-order range
Number_Of_Optimizations = 50; % how many runs of the optimization loop to do (how many initial guesses)
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

%   first, load guess values
    p0 = readmatrix('SSfit2022.csv'); % 2022/23 parameter values (based on HUVECs)
%     p0 = readmatrix('LWCfit2015.csv'); % 2015 parameter values (based on PAECs)

for OptiLoop = 1:Number_Of_Optimizations         
    
    p1 = 10.^((rand(Number_Of_OptimizingParams,1)*5)-6); % initial guesses for all optimizing parameters from 1e-6 to 1e-1
    r1 = 10.^((rand(Number_Of_OptimizingParams,1)*2)-1); % initial guesses for all optimizing parameters with values in x10 - /10 neighborhood

    if InitialGuessType == 1
%       Instead of full random range, here make initial guesses for all
%       optimizing parameters within one order of magnitude of guess above
        p1 = p0(1:15).*r1; % randomize initial guesses
        p1(p1>1e-1) = 1e-1; % upper limit
        p1(p1<1e-6) = 1e-6; % lower limit   
    end
    
    InitGuesses(:,OptiLoop) = p1; 
end

csvwrite('InitialGuessValues.csv',InitGuesses);

%% MAIN LOOP

parfor OptiLoop = 1:Number_Of_Optimizations         
        
    p0 = readmatrix('SSfit2022.csv'); % 2022/23 parameter values (based on HUVECs)
    p1 = InitGuesses(:,OptiLoop); 

    lb=(p1*0)+(1E-6); % lower bound for optimization
    ub=(p1*0)+(1E-1); % upper bound for optimization
    
    % Production rates - assign initial guesses (optimized in nested codes)
    runparams = parameters;
    runparams.kR1prod  = 4.123; %1.73e-4; %2e-3;
    runparams.kR2prod  = 1.058; %1.79e-5; %2e-3;
    runparams.kN1prod  = 4.121; %1.71e-3; %3.3e-3;
    if InitialGuessType == 1
        runparams.kR1prod = p0(16);
        runparams.kR2prod = p0(17);
        runparams.kN1prod = p0(18);        
    end
   
    % Run optimization of unliganded parameters
    options = optimoptions('lsqnonlin','Display','final'); % display output
    optimal = lsqnonlin(@CostFxn_unligandedParams,p1,lb,ub,options, runparams, model, recsurftargets, perturbs); 
    % optimal => optimal values of p1

    % Assign new optimized values of parameters to the parameter set
    runparams.kR1Rab5a     = optimal(1);
    runparams.kR1Rab4at7a  = optimal(2);
    runparams.kR1Rab4a     = optimal(3);
    runparams.kR1Rab4at11a = optimal(4);
    runparams.kR1Rab11a    = optimal(5);

    runparams.kR1N1Rab5a     = optimal(1);
    runparams.kR1N1Rab4at7a  = optimal(2); 
    runparams.kR1N1Rab4a     = optimal(3); 
    runparams.kR1N1Rab4at11a = optimal(4); 
    runparams.kR1N1Rab11a    = optimal(5); 

    runparams.kR2Rab5a     = optimal(6);         
    runparams.kR2Rab4at7a  = optimal(7);     
    runparams.kR2Rab4a     = optimal(8);    
    runparams.kR2Rab4at11a = optimal(9);     
    runparams.kR2Rab11a    = optimal(10);

    runparams.kN1Rab5a     = optimal(11);         
    runparams.kN1Rab4at7a  = optimal(12);     
    runparams.kN1Rab4a     = optimal(13);    
    runparams.kN1Rab4at11a = optimal(14);     
    runparams.kN1Rab11a    = optimal(15); 

    % FINAL RUN WITH OPTIMIZED PARAMETERS

    % after optimizing the trafficking parameters, 
    % optimize the receptor production rates one last time
    p2 = zeros(3,1);
    p2(1)=runparams.kR1prod;
    p2(2)=runparams.kR2prod;
    p2(3)=runparams.kN1prod; 
    initguesses = p2;
    lb=p2*1e-4;     % lower bound
    ub=p2*1e4;     % upper bound
    options = optimoptions('lsqnonlin','Display','off'); % display output
    optimalR = lsqnonlin(@CostFxn_ReceptorProduction,p2,lb,ub,options,recsurftargets,runparams,model);

    runparams.kR1prod = optimalR(1);
    runparams.kR2prod = optimalR(2);
    runparams.kN1prod = optimalR(3);

    OptiLoop
    optimalALL = [optimal;optimalR]
    FinalOptimizations(:,OptiLoop) = optimalALL;

    fileNo = OptiLoop + NumberToBurn;
    parSave(sprintf('outputvals%d.csv',OptiLoop),optimalALL);
    
%     save('InitialGuesses_and_OptimizedValues.mat','InitGuesses','FinalOptimizations');
%     csvwrite('InitialGuessValues.csv',InitGuesses);
%     csvwrite('OptimizedValues.csv',FinalOptimizations);

end

delete(poolobj);
csvwrite('InitialGuessValuesAll.csv',InitGuesses);
csvwrite('OptimizedValuesAll.csv',FinalOptimizations);
