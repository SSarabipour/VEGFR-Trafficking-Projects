function [timepoints, species_out, observables_out] = SimToSteadyState(parameters,runModel);   
% Code running a simulation to steady state
% Typically used for pre-ligand simulation to to get all receptors to steady-state
%
% INPUTS: parameters = parameter set for simulation (gets passed to the model simulation code)
%         runModel   = name of model simulation function code
%         timestp    = spacing of time steps to be returned by the simulation code (not the simulation spacing)
%
% OUTPUTS: timepoints  = simulation time points returned
%          species_out = concentrations of all molecules at each timepoint
%          observables_out = aggregated concentrations for selected metrics

tstart = tic; % start the clock to record time taken to reach steady state
TimeLen  = 24*60*60; % 24 hours, units in seconds; 60*60 = 1 hour
TimeStep = 3*60*60;
Timepoints  = 0:TimeStep:TimeLen; 

%% Initialize simulation - generic initial conditions - simulate day 1
[~, timepoints, species_out, observables_out] = eval(append(runModel,"(Timepoints', [], parameters, 1)"));
daysToStSt = 1;  % counting how many units of time needed to converge to steady state
TimeStart = TimeLen; % starting point for day 2 is the end of day 1

%% Simulate next day, and check convergence to steady state - loop if not converged
diff = 1; 
while (diff>0.0001 & daysToStSt<1000)   % criterion for convergence (and criterion for bailing out of overproduction scenarios)
    % run day n+1
    species_in=species_out(end,:); % end-of-sim concentrations to serve as initial conditions for next sim
    [~, timepoints2, species_out2, observables_out2] = eval(append(runModel,"(Timepoints', species_in, parameters, 1)")); % run day n+1
    % calculate convergence, day n+1 vs day n
    startSp = species_in;
    endSp = species_out2(end,:);
    
    % remove Rab7 elements (because these don't converge!)
    if runModel == "Unligated_VEGFR_model_20230831" % receptor-only model
        Rab7elements = [18:20 26:28 31:32];
    end
    startSp(Rab7elements) = 0.0; 
    endSp(Rab7elements) = 0.0; 
    
    % ignore tiny numbers - Matlab issue
    lowThreshold = (max(species_in)/1e10); 
    startSp(startSp<lowThreshold) = 0.0;
    endSp(endSp<lowThreshold) = 0.0;
        
    % evaluate max difference
    diffSp = species_in*0; % initialize array
    diffSp = (endSp - startSp)./startSp; % fractional change from previous day
    diffSp(isnan(diffSp)) = 0;
    [diff,index] = max(abs(diffSp)); % identify furthest from steady state

    % accumulate the metrics of record 
    timepoints = [timepoints;timepoints2(2:end)+TimeStart];
    species_out = [species_out;species_out2(2:end,:)];
    observables_out = [observables_out;observables_out2(2:end,:)];
    daysToStSt = daysToStSt+1;
    if(mod(daysToStSt,100)==0)
        daysToStSt
        index
    end
    TimeStart = TimeLen*daysToStSt;
end

tend = toc(tstart); % stop the clock to record time taken to reach steady state
% simtoststdetails = table(daysToStSt,tend, diff, index) % output report on steady state to command window



