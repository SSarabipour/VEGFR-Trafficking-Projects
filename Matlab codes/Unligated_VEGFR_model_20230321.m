function [err, timepoints, species_out, observables_out] = Unligated_VEGFR_model_20230321( timepoints, species_init, params, suppress_plot )
% FIX CALLS: function [timepoints, species_out, observables_out] = Unligated_VEGFR_model_20230321(timeLen, timestp, parameters, species_in)
% Generated by BioNetGen based on model built by Sarvenaz Sarabipour PhD and Feilim Mac Gabhann PhD
% Affiliation: Johns Hopkins University, Maryland, United States

%UNLIGATED_VEGFR_MODEL_20230321 Integrate reaction network and plot observables.
%   Integrates the reaction network corresponding to the BioNetGen model
%   'Unligated_VEGFR_model_20230321' and then (optionally) plots the observable trajectories,
%   or species trajectories if no observables are defined. Trajectories are
%   generated using either default or user-defined parameters and initial
%   species values. Integration is performed by the MATLAB stiff solver
%   'ode15s'. UNLIGATED_VEGFR_MODEL_20230321 returns an error value, a vector of timepoints,
%   species trajectories, and observable trajectories.
%   
%   [err, timepoints, species_out, observables_out]
%        = Unligated_VEGFR_model_20230321( timepoints, species_init, parameters, suppress_plot )
%
%   INPUTS:
%   -------
%   species_init    : row vector of 32 initial species populations.
%   timepoints      : column vector of time points returned by integrator.
%   parameters      : row vector of 32 model parameters.
%   suppress_plot   : 0 if a plot is desired (default), 1 if plot is suppressed.
%
%   Note: to specify default value for an input argument, pass the empty array.
%
%   OUTPUTS:
%   --------
%   err             : 0 if the integrator exits without error, non-zero otherwise.
%   timepoints      : a row vector of timepoints returned by the integrator.
%   species_out     : array of species population trajectories
%                        (columns correspond to species, rows correspond to time).
%   observables_out : array of observable trajectories
%                        (columns correspond to observables, rows correspond to time).
%
% convert between params (structure) and parameters (indexed) in the correct order

a = fileread("Paramlist.txt");
pn = strsplit(a);
for i=1:length(pn)
    parameters(i) = params.(pn{i});
end

%% Process input arguments

% define any missing arguments
if ( nargin < 1 )
    timepoints = [];
end

if ( nargin < 2 )
    species_init = [];
end

if ( nargin < 3 )
    parameters = [];
end

if ( nargin < 4 )
    suppress_plot = 0;
end


% initialize outputs (to avoid error msgs if script terminates early
err = 0;
species_out     = [];
observables_out = [];


% setup default parameters, if necessary
if ( isempty(parameters) )
   parameters = [ 1.8, 4.9, 68, 0.008, 0.01, 0.002, 0.01, 0.00166, 0.01, 0.0406440383, 0.0036000689, 0.0025109826, 0.0406440383, 0.00252451198, 0.00261335768, 0.01939407357, 0.00252451198, 0.002163012166, 0.000760959861, 0.019364363943, 0.002163012166, 0.04969202426, 0.01186345288, 0.02087014426, 0.04969202426, 0.00025147487, 0.00031683254, 0.00000551607, 0.00025147487, 3.76E-3, 1.52E-3, 3.74E-3 ];
end
% check that parameters has proper dimensions
if (  size(parameters,1) ~= 1  ||  size(parameters,2) ~= 32  )
    fprintf( 1, 'Error: size of parameter argument is invalid! Correct size = [1 32].\n' );
    err = 1;
    return;
end

% setup default initial values, if necessary
if ( isempty(species_init) )
   species_init = initialize_species( parameters );
end
% check that species_init has proper dimensions
if (  size(species_init,1) ~= 1  ||  size(species_init,2) ~= 32  )
    fprintf( 1, 'Error: size of species_init argument is invalid! Correct size = [1 32].\n' );
    err = 1;
    return;
end

% setup default timepoints, if necessary
if ( isempty(timepoints) )
   timepoints = linspace(0,1800,10000+1)';
end
% check that timepoints has proper dimensions
if (  size(timepoints,1) < 2  ||  size(timepoints,2) ~= 1  )
    fprintf( 1, 'Error: size of timepoints argument is invalid! Correct size = [t 1], t>1.\n' );
    err = 1;
    return;
end

% setup default suppress_plot, if necessary
if ( isempty(suppress_plot) )
   suppress_plot = 0;
end
% check that suppress_plot has proper dimensions
if ( size(suppress_plot,1) ~= 1  ||  size(suppress_plot,2) ~= 1 )
    fprintf( 1, 'Error: suppress_plots argument should be a scalar!\n' );
    err = 1;
    return;
end

% define parameter labels (this is for the user's reference!)
param_labels = { 'VEGFR1_0', 'VEGFR2_0', 'NRP1_0', 'kR1R1on', 'kR1R1off', 'kR2R2on', 'kR2R2off', 'kN1R1on', 'kN1R1off', 'kR1Rab5a', 'kR2Rab5a', 'kN1Rab5a', 'kR1N1Rab5a', 'kR1Rab4a', 'kR2Rab4a', 'kN1Rab4a', 'kR1N1Rab4a', 'kR1Rab4at11a', 'kR2Rab4at11a', 'kN1Rab4at11a', 'kR1N1Rab4at11a', 'kR1Rab11a', 'kR2Rab11a', 'kN1Rab11a', 'kR1N1Rab11a', 'kR1Rab4at7a', 'kR2Rab4at7a', 'kN1Rab4at7a', 'kR1N1Rab4at7a', 'kR1prod', 'kR2prod', 'kN1prod' };



%% Integrate Network Model
 
% calculate expressions
[expressions] = calc_expressions( parameters );

% set ODE integrator options
opts = odeset( 'RelTol',   1e-8,   ...
               'AbsTol',   1e-8,   ...
               'Stats',    'off',  ...
               'BDF',      'off',    ...
               'MaxOrder', 5   );


% define derivative function
rhs_fcn = @(t,y)( calc_species_deriv( t, y, expressions ) );

% simulate model system (stiff integrator)
try 
    [~, species_out] = ode15s( rhs_fcn, timepoints, species_init', opts );
    if(length(timepoints) ~= size(species_out,1))
        exception = MException('ODE15sError:MissingOutput','Not all timepoints output\n');
        throw(exception);
    end
catch
    err = 1;
    fprintf( 1, 'Error: some problem encountered while integrating ODE network!\n' );
    return;
end

% calculate observables
observables_out = zeros( length(timepoints), 40 );
for t = 1 : length(timepoints)
    observables_out(t,:) = calc_observables( species_out(t,:), expressions );
end


%% Plot Output, if desired

if ( ~suppress_plot )
    
    % define plot labels
    observable_labels = { 'R1tot', 'R2tot', 'N1tot', 'R1tots', 'R2tots', 'N1tots', 'R1totRab4a5a', 'R2totRab4a5a', 'N1totRab4a5a', 'R1totRab11a', 'R2totRab11a', 'N1totRab11a', 'R1totRab7a', 'R2totRab7a', 'N1totRab7a', 'R1R1tot', 'R2R2tot', 'R1N1tot', 'R1R1N1tot', 'N1R1R1N1tot', 'R1R1tots', 'R2R2tots', 'R1N1tots', 'R1R1N1tots', 'N1R1R1N1tots', 'R1R1totRab4a5a', 'R2R2totRab4a5a', 'R1N1totRab4a5a', 'R1R1N1totRab4a5a', 'N1R1R1N1totRab4a5a', 'R1R1totRab11a', 'R2R2totRab11a', 'R1N1totRab11a', 'R1R1N1totRab11a', 'N1R1R1N1totRab11a', 'R1R1totRab7a', 'R2R2totRab7a', 'R1N1totRab7a', 'R1R1N1totRab7a', 'N1R1R1N1totRab7a' };

    % construct figure
    plot(timepoints,observables_out);
    title('Unligated_VEGFR_model_Bionetgen observables','fontSize',14,'Interpreter','none');
    axis([0 timepoints(end) 0 inf]);
    legend(observable_labels,'fontSize',10,'Interpreter','none');
    xlabel('time','fontSize',12,'Interpreter','none');
    ylabel('number or concentration','fontSize',12,'Interpreter','none');

end

%~~~~~~~~~~~~~~~~~~~~~%
% END of main script! %
%~~~~~~~~~~~~~~~~~~~~~%

% Define if function to allow nested if statements in user-defined functions
function [val] = if__fun (cond, valT, valF)
% IF__FUN Select between two possible return values depending on the boolean
% variable COND.
    if (cond)
        val = valT;
    else
        val = valF;
    end
end

% initialize species function
function [species_init] = initialize_species( params )

    species_init = zeros(1,32);
    species_init(1) = params(1);
    species_init(2) = params(2);
    species_init(3) = params(3);
    species_init(4:32) = 0;

end


% Calculate expressions
function [ expressions ] = calc_expressions ( parameters )

    expressions = zeros(1,32);
    expressions(1) = parameters(1);
    expressions(2) = parameters(2);
    expressions(3) = parameters(3);
    expressions(4) = parameters(4);
    expressions(5) = parameters(5);
    expressions(6) = parameters(6);
    expressions(7) = parameters(7);
    expressions(8) = parameters(8);
    expressions(9) = parameters(9);
    expressions(10) = parameters(10);
    expressions(11) = parameters(11);
    expressions(12) = parameters(12);
    expressions(13) = parameters(13);
    expressions(14) = parameters(14);
    expressions(15) = parameters(15);
    expressions(16) = parameters(16);
    expressions(17) = parameters(17);
    expressions(18) = parameters(18);
    expressions(19) = parameters(19);
    expressions(20) = parameters(20);
    expressions(21) = parameters(21);
    expressions(22) = parameters(22);
    expressions(23) = parameters(23);
    expressions(24) = parameters(24);
    expressions(25) = parameters(25);
    expressions(26) = parameters(26);
    expressions(27) = parameters(27);
    expressions(28) = parameters(28);
    expressions(29) = parameters(29);
    expressions(30) = parameters(30);
    expressions(31) = parameters(31);
    expressions(32) = parameters(32);
   
end



% Calculate observables
function [ observables ] = calc_observables ( species, expressions )

    observables = zeros(1,40);
    observables(1) = species(1) +2*species(4) +species(6) +species(7) +2*species(10) +2*species(11) +2*species(12) +species(14) +species(15) +species(18) +2*species(21) +2*species(22) +2*species(23) +species(25) +2*species(26) +species(28) +2*species(29) +2*species(30) +2*species(31) +2*species(32);
    observables(2) = species(2) +2*species(5) +species(8) +2*species(13) +species(16) +species(19) +2*species(24) +2*species(27);
    observables(3) = species(3) +species(6) +species(9) +species(10) +2*species(11) +species(14) +species(17) +species(20) +species(21) +2*species(22) +species(25) +species(28) +species(29) +2*species(30) +species(31) +2*species(32);
    observables(4) = species(1) +2*species(4) +species(6) +2*species(10) +2*species(11);
    observables(5) = species(2) +2*species(5);
    observables(6) = species(3) +species(6) +species(10) +2*species(11);
    observables(7) = species(7) +2*species(12) +species(14) +2*species(21) +2*species(22);
    observables(8) = species(8) +2*species(13);
    observables(9) = species(9) +species(14) +species(21) +2*species(22);
    observables(10) = species(15) +2*species(23) +species(25) +2*species(29) +2*species(30);
    observables(11) = species(16) +2*species(24);
    observables(12) = species(17) +species(25) +species(29) +2*species(30);
    observables(13) = species(18) +2*species(26) +species(28) +2*species(31) +2*species(32);
    observables(14) = species(19) +2*species(27);
    observables(15) = species(20) +species(28) +species(31) +2*species(32);
    observables(16) = 2*species(4) +2*species(10) +2*species(11) +2*species(12) +2*species(21) +2*species(22) +2*species(23) +2*species(26) +2*species(29) +2*species(30) +2*species(31) +2*species(32);
    observables(17) = 2*species(5) +2*species(13) +2*species(24) +2*species(27);
    observables(18) = species(6) +2*species(10) +4*species(11) +species(14) +2*species(21) +4*species(22) +species(25) +species(28) +2*species(29) +4*species(30) +2*species(31) +4*species(32);
    observables(19) = 2*species(10) +4*species(11) +2*species(21) +4*species(22) +2*species(29) +4*species(30) +2*species(31) +4*species(32);
    observables(20) = 4*species(11) +4*species(22) +4*species(30) +4*species(32);
    observables(21) = 2*species(4) +2*species(10) +2*species(11);
    observables(22) = 2*species(5);
    observables(23) = species(6) +2*species(10) +4*species(11);
    observables(24) = 2*species(10) +4*species(11);
    observables(25) = 4*species(11);
    observables(26) = 2*species(12) +2*species(21) +2*species(22);
    observables(27) = 2*species(13);
    observables(28) = species(14) +2*species(21) +4*species(22);
    observables(29) = 2*species(21) +4*species(22);
    observables(30) = 4*species(22);
    observables(31) = 2*species(23) +2*species(29) +2*species(30);
    observables(32) = 2*species(24);
    observables(33) = species(25) +2*species(29) +4*species(30);
    observables(34) = 2*species(29) +4*species(30);
    observables(35) = 4*species(30);
    observables(36) = 2*species(26) +2*species(31) +2*species(32);
    observables(37) = 2*species(27);
    observables(38) = species(28) +2*species(31) +4*species(32);
    observables(39) = 2*species(31) +4*species(32);
    observables(40) = 4*species(32);

end


% Calculate ratelaws
function [ ratelaws ] = calc_ratelaws ( species, expressions, observables )

    ratelaws = zeros(1,57);
    ratelaws(1) = 0.5*expressions(4)*species(1)*species(1);
    ratelaws(2) = 0.5*expressions(6)*species(2)*species(2);
    ratelaws(3) = expressions(8)*species(1)*species(3);
    ratelaws(4) = expressions(10)*species(1);
    ratelaws(5) = expressions(11)*species(2);
    ratelaws(6) = expressions(12)*species(3);
    ratelaws(7) = expressions(30);
    ratelaws(8) = expressions(31);
    ratelaws(9) = expressions(32);
    ratelaws(10) = expressions(4)*species(1)*species(6);
    ratelaws(11) = 0.5*expressions(4)*species(6)*species(6);
    ratelaws(12) = expressions(5)*species(4);
    ratelaws(13) = expressions(7)*species(5);
    ratelaws(14) = 2*expressions(8)*species(4)*species(3);
    ratelaws(15) = expressions(9)*species(6);
    ratelaws(16) = expressions(10)*species(4);
    ratelaws(17) = expressions(11)*species(5);
    ratelaws(18) = expressions(13)*species(6);
    ratelaws(19) = expressions(14)*species(7);
    ratelaws(20) = expressions(15)*species(8);
    ratelaws(21) = expressions(16)*species(9);
    ratelaws(22) = expressions(18)*species(7);
    ratelaws(23) = expressions(19)*species(8);
    ratelaws(24) = expressions(20)*species(9);
    ratelaws(25) = expressions(26)*species(7);
    ratelaws(26) = expressions(27)*species(8);
    ratelaws(27) = expressions(28)*species(9);
    ratelaws(28) = expressions(5)*species(10);
    ratelaws(29) = expressions(5)*species(11);
    ratelaws(30) = expressions(8)*species(10)*species(3);
    ratelaws(31) = expressions(9)*species(10);
    ratelaws(32) = 2*expressions(9)*species(11);
    ratelaws(33) = expressions(13)*species(10);
    ratelaws(34) = expressions(13)*species(11);
    ratelaws(35) = expressions(14)*species(12);
    ratelaws(36) = expressions(15)*species(13);
    ratelaws(37) = expressions(17)*species(14);
    ratelaws(38) = expressions(18)*species(12);
    ratelaws(39) = expressions(19)*species(13);
    ratelaws(40) = expressions(21)*species(14);
    ratelaws(41) = expressions(22)*species(15);
    ratelaws(42) = expressions(23)*species(16);
    ratelaws(43) = expressions(24)*species(17);
    ratelaws(44) = expressions(26)*species(12);
    ratelaws(45) = expressions(27)*species(13);
    ratelaws(46) = expressions(29)*species(14);
    ratelaws(47) = expressions(17)*species(21);
    ratelaws(48) = expressions(17)*species(22);
    ratelaws(49) = expressions(21)*species(21);
    ratelaws(50) = expressions(21)*species(22);
    ratelaws(51) = expressions(22)*species(23);
    ratelaws(52) = expressions(23)*species(24);
    ratelaws(53) = expressions(25)*species(25);
    ratelaws(54) = expressions(29)*species(21);
    ratelaws(55) = expressions(29)*species(22);
    ratelaws(56) = expressions(25)*species(29);
    ratelaws(57) = expressions(25)*species(30);

end

% Calculate species derivatives
function [ Dspecies ] = calc_species_deriv ( time, species, expressions )
    
    % initialize derivative vector
    Dspecies = zeros(32,1);
    
    % update observables
    [ observables ] = calc_observables( species, expressions );
    
    % update ratelaws
    [ ratelaws ] = calc_ratelaws( species, expressions, observables );
                        
    % calculate derivatives
    Dspecies(1) = -2.0*ratelaws(1) -ratelaws(3) -ratelaws(4) +ratelaws(7) -ratelaws(10) +2.0*ratelaws(12) +ratelaws(15) +ratelaws(19) +ratelaws(28) +ratelaws(41);
    Dspecies(2) = -2.0*ratelaws(2) -ratelaws(5) +ratelaws(8) +2.0*ratelaws(13) +ratelaws(20) +ratelaws(42);
    Dspecies(3) = -ratelaws(3) -ratelaws(6) +ratelaws(9) -ratelaws(14) +ratelaws(15) +ratelaws(21) -ratelaws(30) +ratelaws(31) +ratelaws(32) +ratelaws(43);
    Dspecies(4) = ratelaws(1) -ratelaws(12) -ratelaws(14) -ratelaws(16) +ratelaws(31) +ratelaws(35) +ratelaws(51);
    Dspecies(5) = ratelaws(2) -ratelaws(13) -ratelaws(17) +ratelaws(36) +ratelaws(52);
    Dspecies(6) = ratelaws(3) -ratelaws(10) -2.0*ratelaws(11) -ratelaws(15) -ratelaws(18) +ratelaws(28) +2.0*ratelaws(29) +ratelaws(37) +ratelaws(53);
    Dspecies(7) = ratelaws(4) -ratelaws(19) -ratelaws(22) -ratelaws(25);
    Dspecies(8) = ratelaws(5) -ratelaws(20) -ratelaws(23) -ratelaws(26);
    Dspecies(9) = ratelaws(6) -ratelaws(21) -ratelaws(24) -ratelaws(27);
    Dspecies(10) = ratelaws(10) +ratelaws(14) -ratelaws(28) -ratelaws(30) -ratelaws(31) +ratelaws(32) -ratelaws(33) +ratelaws(47) +ratelaws(56);
    Dspecies(11) = ratelaws(11) -ratelaws(29) +ratelaws(30) -ratelaws(32) -ratelaws(34) +ratelaws(48) +ratelaws(57);
    Dspecies(12) = ratelaws(16) -ratelaws(35) -ratelaws(38) -ratelaws(44);
    Dspecies(13) = ratelaws(17) -ratelaws(36) -ratelaws(39) -ratelaws(45);
    Dspecies(14) = ratelaws(18) -ratelaws(37) -ratelaws(40) -ratelaws(46);
    Dspecies(15) = ratelaws(22) -ratelaws(41);
    Dspecies(16) = ratelaws(23) -ratelaws(42);
    Dspecies(17) = ratelaws(24) -ratelaws(43);
    Dspecies(18) = ratelaws(25);
    Dspecies(19) = ratelaws(26);
    Dspecies(20) = ratelaws(27);
    Dspecies(21) = ratelaws(33) -ratelaws(47) -ratelaws(49) -ratelaws(54);
    Dspecies(22) = ratelaws(34) -ratelaws(48) -ratelaws(50) -ratelaws(55);
    Dspecies(23) = ratelaws(38) -ratelaws(51);
    Dspecies(24) = ratelaws(39) -ratelaws(52);
    Dspecies(25) = ratelaws(40) -ratelaws(53);
    Dspecies(26) = ratelaws(44);
    Dspecies(27) = ratelaws(45);
    Dspecies(28) = ratelaws(46);
    Dspecies(29) = ratelaws(49) -ratelaws(56);
    Dspecies(30) = ratelaws(50) -ratelaws(57);
    Dspecies(31) = ratelaws(54);
    Dspecies(32) = ratelaws(55);

end


end