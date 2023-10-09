function CostOut = CostFxn_unligandedParams(p1,parameters,model,recsurftargets,perturbs)
%
% This cost function evaluates the experimental fit at the current guess (p1)
% The simulations of experiments are actually in the called function UnligandedHUVECExperimentsAll
% Much of this cost function is setting up the basis for those simulations.
%

tstartcost = tic;

suppressligands = 1;

%% Incorporate current guesses into the current parameter set
parameters.kR1Rab5a     = p1(1);    %This is the R1 internalization rate from surface to Rab5a that we are trying to optimize
parameters.kR1Rab4at7a  = p1(2);    %This is the R1 transfer from Rab4a to Rab7a rate that we are trying to optimize
parameters.kR1Rab4a     = p1(3);    %This is the R1 recycling rate from Rab4a to surface rate that we are trying to optimize
parameters.kR1Rab4at11a = p1(4);    %This is the R1 transfer from Rab4a to Rab11a rate that we are trying to optimize
parameters.kR1Rab11a    = p1(5);    %This is the R1 recycling rate from Rab11a to surface rate that we are trying to optimize

parameters.kR1N1Rab5a     = p1(1);
parameters.kR1N1Rab4at7a  = p1(2); 
parameters.kR1N1Rab4a     = p1(3); 
parameters.kR1N1Rab4at11a = p1(4); 
parameters.kR1N1Rab11a    = p1(5); 

parameters.kR2Rab5a     = p1(6);    %This is the R2 internalization rate from surface to Rab5a that we are trying to optimize
parameters.kR2Rab4at7a  = p1(7);    %This is the R2 transfer from Rab4a to Rab7a rate that we are trying to optimize
parameters.kR2Rab4a     = p1(8);    %This is the R2 recycling rate from Rab4a to surface rate that we are trying to optimize
parameters.kR2Rab4at11a = p1(9);    %This is the R2 transfer from Rab4a to Rab11a rate that we are trying to optimize
parameters.kR2Rab11a    = p1(10);   %This is the R2 recycling rate from Rab11a to surface rate that we are trying to optimize

parameters.kN1Rab5a     = p1(11);    %This is the N1 internalization rate from surface to Rab5a that we are trying to optimize
parameters.kN1Rab4at7a  = p1(12);    %This is the N1 transfer from Rab4a to Rab7a rate that we are trying to optimize
parameters.kN1Rab4a     = p1(13);    %This is the N1 recycling rate from Rab4a to surface rate that we are trying to optimize
parameters.kN1Rab4at11a = p1(14);    %This is the N1 transfer from Rab4a to Rab11a rate that we are trying to optimize
parameters.kN1Rab11a    = p1(15);    %This is the N1 recycling rate from Rab11a to surface rate that we are trying to optimize


%% INNER LOOP - optimizing for receptor production rates

% The Rab4, Rab11 baselines will be modified by knockdown; 
%   note their current values
kR1Rab4a_baseline  = parameters.kR1Rab4a;
kR1Rab11a_baseline = parameters.kR1Rab11a;
        
tstartR = tic; % calculate length of simulation time

% initial guesses for receptor production optimization
%   - if available use previous optimization as starting point
p2(1)=parameters.kR1prod;
p2(2)=parameters.kR2prod;
p2(3)=parameters.kN1prod;
initguesses = p2;
lb=p2*1e-4;     % lower bound
ub=p2*1e4;     % upper bound
options = optimoptions('lsqnonlin','Display','off'); % display output
optimalR = lsqnonlin(@CostFxn_ReceptorProduction,p2,lb,ub,options,recsurftargets,parameters,model);

tend = toc(tstartR);

% output the optimized receptor production rates and other info
optRdetails = table(tend, initguesses, optimalR);
% update the parameter set for this (outer) parameter loop
parameters.kR1prod = optimalR(1);
parameters.kR2prod = optimalR(2);
parameters.kN1prod = optimalR(3);


%% INITIAL STEADY STATE USING THE OPTIMIZED RECEPTOR PRODUCTION RATES
[timepoints, species_out, observables_out] = SimToSteadyState(parameters,model);
ConcOut=CalcOutputs(observables_out,suppressligands);
ConcOutNames = fieldnames(ConcOut);
TotPanels = length(ConcOutNames);
PanelRows = ceil(length(ConcOutNames)/3);


%% OUTER LOOP - UNLIGANDED TRAFFICKING PARAMETERS

speciesInit = species_out(end,:); % we just ran everything to steady state, so use the outputs as new inputs
baselines = CalcBaselines(ConcOut);
% pout = CalcPout(parameters,model,timestp_stst,suppressligands); % basically a list of key values. BEWARE - runs steady state, slows everything down.


%% RUN Unliganded EXPERIMENTS
createFigs = 0;
[R1_exp,R1_sim,R2_exp,R2_sim,N1_exp,N1_sim] = UnligandedHUVECExperimentsAll(parameters,speciesInit,baselines,perturbs,recsurftargets,model,createFigs);

%% COST FUNCTION CALCULATION
% each experiment returns a simulated and experimental value, for each of
% VEGFR1, VEGFR2, and NRP1. Many of these are zeroes for both (as in, no
% data for that receptor for that experiment)

% The error (CostOut) is the difference between the simulated and
% experimental values

% Cost Function type (a) - absolute difference
% CostOut = [R1_sim - R1_exp; R2_sim - R2_exp; N1_sim - N1_exp]; %3*length(R1_theor) elements

% Cost Function type (b) - normalization
                CostOut = [(R1_sim - R1_exp)./R1_exp; (R2_sim - R2_exp)./R2_exp; (N1_sim - N1_exp)./N1_exp]; %3*length(R1_theor) elements
                CostOut (isnan(CostOut)) = 0; 

% Cost Function type (c) fold change - max increase or decrease
%                 CostOuta = [R1_sim./R1_exp; R2_sim./R2_exp; N1_sim./N1_exp]; %3*length(R1_theor) elements
%                 CostOuta(isnan(CostOuta)) = 0;
%                 CostOutb = [R1_exp./R1_sim; R2_exp./R2_sim; N1_exp./N1_sim]; %3*length(R1_theor) elements
%                 CostOutb(isnan(CostOutb)) = 0;
%                 CostOut = max(CostOuta,CostOutb);

%% COST ADJUSTMENTS (if desired)

% % Emphasize the biotinylayion surface-internal simulations
% j=9;
% CostOut(j) = CostOut(j)*10; % give these higher weight to compensate for fewer measurements
% CostOut(j+length(R1_sim)) = CostOut(j+length(R1_sim))*10; 
% CostOut(j+length(R1_sim)+length(R2_sim)) = CostOut(j+length(R1_sim)+length(R2_sim))*10; 
% 
% % Emphasize the siRNA simulations
% j=6:8;
% CostOut(j) = CostOut(j)*3; % give these higher weight to compensate for fewer measurements
% CostOut(j+length(R1_sim)) = CostOut(j+length(R1_sim))*3; 
% CostOut(j+length(R1_sim)+length(R2_sim)) = CostOut(j+length(R1_sim)+length(R2_sim))*3; 

% Error suppression (for debugging purposes)
% Choose which experiments to suppress in cost function:
% Weighting these experiments to zero
% EXAMPLE: suppress1 = [1:5 9]; suppress2 = suppress1+length(R1_sim); suppress = [suppress1 suppress2]; CostOut(suppress) = 0;     

tendcost = toc(tstartcost);

end
