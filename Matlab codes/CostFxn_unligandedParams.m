function CostOut = CostFxn_unligandedParams(p1,parameters,model,recsurftargets,perturbs)
%
% This cost function evaluates the experimental fit at the current guess (p1)
% The simulations of experiments are actually in the called function UnligandedHUVECExperimentsAll
% Much of this cost function is setting up the basis for those simulations.
% Rab7a is the degradation compartment

global CurrGuessFig
global ig

save('current_params.mat','p1')
figure(CurrGuessFig);
    scatter([1:15],ig,'k.');
    hold on;
    scatter([1:15],p1,'bd');
    set(gca,'YScale','log');
    ylim([1e-6 1e-1]);
    drawnow;

tstartcost = tic;

suppressligands = 1;

%% Global Variables (see driver for details)
global r1loop        % rloop is the inner, receptor production rate loop
global ploop         % ploop is the outer, trafficking parameters loop
ploop = ploop + 1;
r1loop = 0;
loops = [ploop r1loop] % output the current loops
        
global currentR1prodguess
global currentR2prodguess
global currentN1prodguess

global h1
global CheckStStateAfterRProdOpti


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

% The Rab4a, Rab11a baselines will be modified by knockdown; note their
% current values
kR1Rab4a_baseline  = parameters.kR1Rab4a;
kR1Rab11a_baseline = parameters.kR1Rab11a;
        
tstartR = tic; % calculate length of simulation time

% use immediate previous values from global var
% initial guesses for receptor production optimization - use previous optimization as starting point
p2(1)=currentR1prodguess;
p2(2)=currentR2prodguess;
p2(3)=currentN1prodguess;
initguesses = p2;
lb=p2*.001;     % lower bound
ub=p2*1000;     % upper bound
options = optimoptions('lsqnonlin','Display','iter'); % display output
optimalR = lsqnonlin(@CostFxn_ReceptorProduction,p2,lb,ub,options,recsurftargets,parameters,model);

tend = toc(tstartR);

% output the optimized receptor production rates and other info
optRdetails = table(tend, initguesses, optimalR)
% update the global variable
currentR1prodguess = optimalR(1);
currentR2prodguess = optimalR(2);
currentN1prodguess = optimalR(3);
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

figure(CheckStStateAfterRProdOpti);
for m=1:TotPanels
    eval(append("ax",num2str(m),"=subplot(",num2str(PanelRows),",3,",num2str(m),");"));
    plot(timepoints/3600,ConcOut.(ConcOutNames{m}));
    hold on;
	title(ConcOutNames{m}, 'Interpreter','none');
end


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

% The error (CostOut) is the different between the simulated and
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

% Error suppression (for debugging purposes only)
% Choose which experiments to suppress in cost function:
% Weighting these experiments (Synthetic data) to zero
% EXAMPLE: suppress = [18:25 26:30 31:36]; CostOut(suppress) = 0;
% OR: suppress1 = [19:21]+length(R1_theor); % suppr R2; suppress2 = [28:33]; suppress = [suppress1 suppress2]; CostOut(suppress) = 0; 
     
%% Calculate and visualize goodness-of-fit
ratiosR1 = [R1_sim./R1_exp];
ratiosR1(isnan(ratiosR1)) = 0;
ratiosR2 = [R2_sim./R2_exp];
ratiosR2(isnan(ratiosR2)) = 0;
ratiosN1 = [N1_sim./N1_exp];
ratiosN1(isnan(ratiosN1)) = 0;

relR1 = R1_sim/R1_sim(1);
relR2 = R2_sim/R2_sim(1);
relN1 = N1_sim/N1_sim(1);

% colorize the bars by experiments

BarColorList = zeros(length(R1_sim),3);
for j = 1:5
    BarColorList(j,:) = [1 0 0]; % KQ1 (CHX expts)
end
for j = 6:8
    BarColorList(j,:) = [0 1 0]; % KQ2 (siRNA vs Rab expts)
end
for j = 9:9
    BarColorList(j,:) = [0 0 0]; % Surface Ratios
end      

figure (h1);

ax1=subplot(3,3,1);
b1=bar(ratiosR1);
b1.FaceColor = 'flat';
b1.CData = BarColorList;
title('err ratio - R1');

ax2=subplot(3,3,2);
b2=bar(ratiosR2);
b2.FaceColor = 'flat';
b2.CData = BarColorList;
title('err ratio - R2');

ax3=subplot(3,3,3);
b3=bar(ratiosN1);
b3.FaceColor = 'flat';
b3.CData = BarColorList;
title('err ratio - N1');

ax4=subplot(3,3,4);
b4=bar(relR1);
b4.FaceColor = 'flat';
b4.CData = BarColorList;
title('rel - R1');

ax5=subplot(3,3,5);
b5=bar(relR2);
b5.FaceColor = 'flat';
b5.CData = BarColorList;
title('rel - R2');

ax6=subplot(3,3,6);
b6=bar(relN1);
b6.FaceColor = 'flat';
b6.CData = BarColorList;
title('rel - R1');

ax7=subplot(3,3,7);
b7=bar(CostOut(1:length(R1_sim)));
b7.FaceColor = 'flat';
b7.CData = BarColorList;
title('cost - R1');

ax8=subplot(3,3,8);
b8=bar(CostOut( (length(R1_sim)+1) : (length(R1_sim)+length(R2_sim)) ));
b8.FaceColor = 'flat';
b8.CData = BarColorList;
title('cost - R2');

ax9=subplot(3,3,9);
b9=bar(CostOut( (length(R1_sim)+length(R2_sim)+1) : (length(R1_sim)+length(R2_sim)+length(N1_sim)) ));
b9.FaceColor = 'flat';
b9.CData = BarColorList;
title('cost - N1');

% drawnow;

tendcost = toc(tstartcost)

end
