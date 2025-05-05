close all;
clear all;

model = "All_Lig_Traff_model_Bionetgen_20231016"; %"TraffPhosph"
parameters = base_parameters(model);

% HIGHER VEGFR2 INTERNALIZATION
parameters.kVegfR2Rab5a = parameters.kR2Rab5a * 3;
parameters.kVegfR2N1Rab5a = parameters.kVegfR2Rab5a;

% ngml = [0 0.01 0.025 0.05 0.1 0.25 0.5 1 2.5 5 10 25 50 100 150];
ngml = 50; % standard ligand initial concentration
ligandforthisrun = "v165a"; % alternatives "v121a"; "plgf1"; "plgf2";

paramadj = [0.2 0.5 1.0 2.0 5.0]; % for adjusting trafficking parameters
paramadjkintR2 = [1.0/3.0 2.0/3.0 3.0/3.0 5.0/3.0 8.0/3.0];

timepts = 0:60:(240*60); %[0 15 30 60 60 240];

ReFitProductionRates = 0 ; % 1 for yes, 0 for no
R1surf_target = 1800;
R2surf_target = 4900; 
N1surf_target = 68000;

visibility = 'off';
longfigs = 1; 

% runCase = "onlyR1";
% runCase = "onlyR2";
% runCase = "onlyN1";
% runCase = "R1&R2";
% runCase = "R1&N1";
% runCase = "N1&R2";
runCase = "All3R";

% DimType = "121"; % assuming 1:1 binding of ligands to receptors
% DimType = "LID"; % assuming no pre-dimerization of receptors
% DimType = "DPD"; % assuming only pre-dimerization of receptors and no ligand-induced dimerization
                   % Note: DPD mode is not quite right (because it still permits ligands binding to R monomers without productive further binding)
DimType = "All"; % all dimerization modes

ScatchardMode = "off"; % when "on" assumes no production, no internalization, only binding (4dC)

%% Surface area of cellular compartments
SA_surface = 1000; % (um2/cell)
SA_rab4    = 950; % (um2/cell)
SA_rab11   = 325; % (um2/cell)

vol_surface = 1e7; % (fL/cell)
vol_rab4    = 11.25; % (fL/cell)
vol_rab11   = 3.75; % (fL/cell)


switch runCase
    case "onlyR1"
        parameters.VEGFR1_0 = R1surf_target; 
        parameters.VEGFR2_0 = 0; 
        parameters.NRP1_0   = 0; 
        parameters.kR2prod  = 0; 
        parameters.kN1prod  = 0; 
        R2surf_target = 0; 
        N1surf_target = 0;
    case "onlyR2"
        parameters.VEGFR1_0 = 0; 
        parameters.VEGFR2_0 = R2surf_target; 
        parameters.NRP1_0   = 0; 
        parameters.kR1prod  = 0; 
        parameters.kN1prod  = 0; 
        R1surf_target = 0;
        N1surf_target = 0;
    case "onlyN1"
        parameters.VEGFR1_0 = 0; 
        parameters.VEGFR2_0 = 0; 
        parameters.NRP1_0   = N1surf_target; 
        parameters.kR1prod  = 0; 
        parameters.kR2prod  = 0; 
        R1surf_target = 0;
        R2surf_target = 0; 
    case "R1&R2"
        parameters.VEGFR1_0 = R1surf_target; 
        parameters.VEGFR2_0 = R2surf_target; 
        parameters.NRP1_0   = 0; 
        parameters.kN1prod  = 0; 
        N1surf_target = 0;
    case "R1&N1"
        parameters.VEGFR1_0 = R1surf_target; 
        parameters.VEGFR2_0 = 0; 
        parameters.NRP1_0   = N1surf_target; 
        parameters.kR2prod  = 0; 
        R2surf_target = 0; 
    case "N1&R2"
        parameters.VEGFR1_0 = 0; 
        parameters.VEGFR2_0 = R2surf_target; 
        parameters.NRP1_0   = N1surf_target; 
        parameters.kR1prod  = 0; 
        R1surf_target = 0;
    case "All3R"
        parameters.VEGFR1_0 = R1surf_target; 
        parameters.VEGFR2_0 = R2surf_target; 
        parameters.NRP1_0   = N1surf_target; 
end

if DimType == "121"
    parameters.VEGFR1_0 = parameters.VEGFR1_0/2; 
    parameters.VEGFR2_0 = parameters.VEGFR2_0/2;  
    parameters.NRP1_0   = parameters.NRP1_0/2; 
end

switch ScatchardMode
    case "on"
        % no production or internalization (receptors)
        parameters.kN1prod  = 0; 
        parameters.kR2prod  = 0; 
        parameters.kR1prod  = 0; 
        
        parameters.kR2Rab5a = 0; 
        parameters.kR1Rab5a = 0; 
        parameters.kN1Rab5a = 0; 
        parameters.kR2N1Rab5a = 0; 
        parameters.kR1N1Rab5a = 0; 

        % no internalization (ligand-receptor complexes)
        parameters.kVegfR1Rab5a = parameters.kR1Rab5a;
        parameters.kVegfR1N1Rab5a = parameters.kVegfR1Rab5a;
        parameters.kPlgfR1Rab5a = parameters.kR1Rab5a;
        parameters.kPlgfR1N1Rab5a = parameters.kPlgfR1Rab5a;
        parameters.kVegfR2Rab5a = parameters.kR2Rab5a * 3;
        parameters.kVegfR2N1Rab5a = parameters.kVegfR2Rab5a;

        parameters.kVegfN1Rab5a = parameters.kN1Rab5a;
        parameters.kPlgfN1Rab5a = parameters.kN1Rab5a;

        parameters.kVegfRab5a = 0;
        parameters.kPlgfRab5a = 0;
end

    
switch DimType
    case "LID"
        % No receptor coupling
        kR1R1on = 0; % units = (receptors/um2)-1.s-1
        parameters.kR1R1on_surf  = kR1R1on/SA_surface;
        parameters.kR1R1on_rab4  = kR1R1on/SA_rab4;
        parameters.kR1R1on_rab11 = kR1R1on/SA_rab11;
        parameters.kR1R1off = 0.01;
        kR2R2on = 0; % units = (receptors/um2)-1.s-1
        parameters.kR2R2on_surf  = kR2R2on/SA_surface;
        parameters.kR2R2on_rab4  = kR2R2on/SA_rab4;
        parameters.kR2R2on_rab11 = kR2R2on/SA_rab11;
        parameters.kR2R2off = 0.01;        
        parameters.kN1R1on_surf = parameters.kR1R1on_surf;
        parameters.kN1R1on_rab4 = parameters.kR1R1on_rab4;
        parameters.kN1R1on_rab11 = parameters.kR1R1on_rab11;
        parameters.kN1R1off = parameters.kR1R1off;
        parameters.kplgfR1_N1on_surf = parameters.kN1R1on_surf;
        parameters.kplgfR1_N1on_rab4 = parameters.kR1R1on_rab4;
        parameters.kplgfR1_N1on_rab11 = parameters.kR1R1on_rab11;
        parameters.kplgfR1_N1off = parameters.kN1R1off;
        parameters.kvegfAR1_N1on_surf = parameters.kN1R1on_surf;
        parameters.kvegfAR1_N1on_rab4 = parameters.kR1R1on_rab4;
        parameters.kvegfAR1_N1on_rab11 = parameters.kR1R1on_rab11;
        parameters.kvegfAR1_N1off = parameters.kN1R1off;

        % remove internal coupling
        parameters.kdeltaR1R1_V = 0;
        parameters.kdeltaR1R1_P = 0;
        parameters.kdeltaR2R2 = 0;
        parameters.kdeltaVEGFA_R = 0;
        parameters.kdeltaPLGF_R = 0;

        % adjustment to R2VR2, R1VR1 coupling rates
        adjcouplrates = 1.0;
        parameters.kR2vegfA_R2on_surf = parameters.kR2vegfA_R2on_surf * adjcouplrates;
        parameters.kR2vegfA_R2on_rab4 = parameters.kR2vegfA_R2on_rab4 * adjcouplrates;
        parameters.kR2vegfA_R2on_rab11 = parameters.kR2vegfA_R2on_rab11 * adjcouplrates;

        parameters.kR1vegfA_R1on_surf = parameters.kR1vegfA_R1on_surf * adjcouplrates;
        parameters.kR1vegfA_R1on_rab4 = parameters.kR1vegfA_R1on_rab4 * adjcouplrates;
        parameters.kR1vegfA_R1on_rab11 = parameters.kR1vegfA_R1on_rab11 * adjcouplrates;
        parameters.kR1plgf_R1on_surf = parameters.kR1plgf_R1on_surf * adjcouplrates;
        parameters.kR1plgf_R1on_rab4 = parameters.kR1plgf_R1on_rab4 * adjcouplrates;
        parameters.kR1plgf_R1on_rab11 = parameters.kR1plgf_R1on_rab11 * adjcouplrates;

        parameters.kN1vegfA_N1on_surf = parameters.kN1vegfA_N1on_surf * adjcouplrates;
        parameters.kN1vegfA_N1on_rab4 = parameters.kN1vegfA_N1on_rab4 * adjcouplrates;
        parameters.kN1vegfA_N1on_rab11 = parameters.kN1vegfA_N1on_rab11 * adjcouplrates;
        parameters.kN1plgf2_N1on_surf = parameters.kN1plgf2_N1on_surf * adjcouplrates;
        parameters.kN1plgf2_N1on_rab4 = parameters.kN1plgf2_N1on_rab4 * adjcouplrates;
        parameters.kN1plgf2_N1on_rab11 = parameters.kN1plgf2_N1on_rab11 * adjcouplrates;
        
    case "DPD"

       % Set RVR coupling to zero
        parameters.kR1vegfA_R1on_surf = 0;
        parameters.kR1vegfA_R1on_rab4 = parameters.kR1vegfA_R1on_surf * SA_surface/SA_rab4;
        parameters.kR1vegfA_R1on_rab11 = parameters.kR1vegfA_R1on_surf * SA_surface/SA_rab11;
        parameters.kR1plgf_R1on_surf = 0;
        parameters.kR1plgf_R1on_rab4 = parameters.kR1plgf_R1on_surf * SA_surface/SA_rab4;
        parameters.kR1plgf_R1on_rab11 = parameters.kR1plgf_R1on_surf * SA_surface/SA_rab11;
        parameters.kR2vegfA_R2on_surf = 0;
        parameters.kR2vegfA_R2on_rab4 = parameters.kR2vegfA_R2on_surf * SA_surface/SA_rab4;
        parameters.kR2vegfA_R2on_rab11 = parameters.kR2vegfA_R2on_surf * SA_surface/SA_rab11;
        parameters.kN1vegfA_N1on_surf = 0;
        parameters.kN1vegfA_N1on_rab4 = parameters.kN1vegfA_N1on_surf * SA_surface/SA_rab4;
        parameters.kN1vegfA_N1on_rab11 = parameters.kN1vegfA_N1on_surf * SA_surface/SA_rab11;
        parameters.kN1vegfA_N1off = parameters.kvegfA_N1off;
        parameters.kN1plgf2_N1on_surf = 0;
        parameters.kN1plgf2_N1on_rab4 = parameters.kN1plgf2_N1on_surf * SA_surface/SA_rab4;
        parameters.kN1plgf2_N1on_rab11 = parameters.kN1plgf2_N1on_surf * SA_surface/SA_rab11;
        parameters.kN1plgf2_N1off = parameters.kplgf_N1off;

        % remove internal coupling
%         parameters.kdeltaR1R1_V = 0;
%         parameters.kdeltaR1R1_P = 0;
%         parameters.kdeltaR2R2 = 0;
%         parameters.kdeltaVEGFA_R = 0;
%         parameters.kdeltaPLGF_R = 0;

    case "All"

        % adjustment to R2VR2, R1VR1 coupling rates
        adjcouplrates = 1.0;
        parameters.kR2vegfA_R2on_surf = parameters.kR2vegfA_R2on_surf * adjcouplrates;
        parameters.kR2vegfA_R2on_rab4 = parameters.kR2vegfA_R2on_rab4 * adjcouplrates;
        parameters.kR2vegfA_R2on_rab11 = parameters.kR2vegfA_R2on_rab11 * adjcouplrates;

        parameters.kR1vegfA_R1on_surf = parameters.kR1vegfA_R1on_surf * adjcouplrates;
        parameters.kR1vegfA_R1on_rab4 = parameters.kR1vegfA_R1on_rab4 * adjcouplrates;
        parameters.kR1vegfA_R1on_rab11 = parameters.kR1vegfA_R1on_rab11 * adjcouplrates;
        parameters.kR1plgf_R1on_surf = parameters.kR1plgf_R1on_surf * adjcouplrates;
        parameters.kR1plgf_R1on_rab4 = parameters.kR1plgf_R1on_rab4 * adjcouplrates;
        parameters.kR1plgf_R1on_rab11 = parameters.kR1plgf_R1on_rab11 * adjcouplrates;

        parameters.kN1vegfA_N1on_surf = parameters.kN1vegfA_N1on_surf * adjcouplrates;
        parameters.kN1vegfA_N1on_rab4 = parameters.kN1vegfA_N1on_rab4 * adjcouplrates;
        parameters.kN1vegfA_N1on_rab11 = parameters.kN1vegfA_N1on_rab11 * adjcouplrates;
        parameters.kN1plgf2_N1on_surf = parameters.kN1plgf2_N1on_surf * adjcouplrates;
        parameters.kN1plgf2_N1on_rab4 = parameters.kN1plgf2_N1on_rab4 * adjcouplrates;
        parameters.kN1plgf2_N1on_rab11 = parameters.kN1plgf2_N1on_rab11 * adjcouplrates;

    case "121"
                
        % No receptor coupling
        kR1R1on = 0; % units = (receptors/um2)-1.s-1
        parameters.kR1R1on_surf  = kR1R1on/SA_surface;
        parameters.kR1R1on_rab4  = kR1R1on/SA_rab4;
        parameters.kR1R1on_rab11 = kR1R1on/SA_rab11;
        parameters.kR1R1off = 0.01;
        kR2R2on = 0; % units = (receptors/um2)-1.s-1
        parameters.kR2R2on_surf  = kR2R2on/SA_surface;
        parameters.kR2R2on_rab4  = kR2R2on/SA_rab4;
        parameters.kR2R2on_rab11 = kR2R2on/SA_rab11;
        parameters.kR2R2off = 0.01;        
        parameters.kN1R1on_surf = parameters.kR1R1on_surf;
        parameters.kN1R1on_rab4 = parameters.kR1R1on_rab4;
        parameters.kN1R1on_rab11 = parameters.kR1R1on_rab11;
        parameters.kN1R1off = parameters.kR1R1off;
        parameters.kplgfR1_N1on_surf = parameters.kN1R1on_surf;
        parameters.kplgfR1_N1on_rab4 = parameters.kR1R1on_rab4;
        parameters.kplgfR1_N1on_rab11 = parameters.kR1R1on_rab11;
        parameters.kplgfR1_N1off = parameters.kN1R1off;
        parameters.kvegfAR1_N1on_surf = parameters.kN1R1on_surf;
        parameters.kvegfAR1_N1on_rab4 = parameters.kR1R1on_rab4;
        parameters.kvegfAR1_N1on_rab11 = parameters.kR1R1on_rab11;
        parameters.kvegfAR1_N1off = parameters.kN1R1off;

        % Replace L-R binding rates with 1:1 binding rates
        parameters.kvegfA_R1on_surf = parameters.kvegfA_R1on_surf * 2;
        parameters.kvegfA_R1on_rab4 = parameters.kvegfA_R1on_surf * vol_surface/vol_rab4;
        parameters.kvegfA_R1on_rab11 = parameters.kvegfA_R1on_surf * vol_surface/vol_rab11;
        parameters.kvegfA_R1off = 0.001;

        parameters.kN1R1_vegfAon_surf = parameters.kvegfA_R1on_surf;
        parameters.kN1R1_vegfAon_rab4 = parameters.kvegfA_R1on_rab4;
        parameters.kN1R1_vegfAon_rab11 = parameters.kvegfA_R1on_rab11;
        parameters.kN1R1_vegfAoff = parameters.kvegfA_R1off;

        parameters.kplgf_R1on_surf = parameters.kplgf_R1on_surf * 2;
        parameters.kplgf_R1on_rab4 = parameters.kplgf_R1on_surf * vol_surface/vol_rab4;
        parameters.kplgf_R1on_rab11 = parameters.kplgf_R1on_surf * vol_surface/vol_rab11;
        parameters.kplgf_R1off = 3.5e-4;

        parameters.kN1R1_plgfon_surf = parameters.kplgf_R1on_surf;
        parameters.kN1R1_plgfon_rab4 = parameters.kplgf_R1on_rab4;
        parameters.kN1R1_plgfon_rab11 = parameters.kplgf_R1on_rab11;
        parameters.kN1R1_plgfoff = parameters.kplgf_R1off;

        parameters.kvegfA_R2on_surf = parameters.kvegfA_R2on_surf * 2;
        parameters.kvegfA_R2on_rab4 = parameters.kvegfA_R2on_surf * vol_surface/vol_rab4;
        parameters.kvegfA_R2on_rab11 = parameters.kvegfA_R2on_surf * vol_surface/vol_rab11;
        parameters.kvegfA_R2off = 0.001;

        parameters.kvegfA_N1on_surf = parameters.kvegfA_N1on_surf * 2;
        parameters.kvegfA_N1on_rab4 = parameters.kvegfA_N1on_surf * vol_surface/vol_rab4;
        parameters.kvegfA_N1on_rab11 = parameters.kvegfA_N1on_surf * vol_surface/vol_rab11;
        parameters.kvegfA_N1off = 6e-4;

        parameters.kplgf_N1on_surf = parameters.kplgf_N1on_surf * 2;
        parameters.kplgf_N1on_rab4 = parameters.kplgf_N1on_surf * vol_surface/vol_rab4;
        parameters.kplgf_N1on_rab11 = parameters.kplgf_N1on_surf * vol_surface/vol_rab11;
        parameters.kplgf_N1off = 1e-3;

        % Set RVR coupling to zero
        parameters.kR1vegfA_R1on_surf = 0;
        parameters.kR1vegfA_R1on_rab4 = parameters.kR1vegfA_R1on_surf * SA_surface/SA_rab4;
        parameters.kR1vegfA_R1on_rab11 = parameters.kR1vegfA_R1on_surf * SA_surface/SA_rab11;
        parameters.kR1plgf_R1on_surf = 0;
        parameters.kR1plgf_R1on_rab4 = parameters.kR1plgf_R1on_surf * SA_surface/SA_rab4;
        parameters.kR1plgf_R1on_rab11 = parameters.kR1plgf_R1on_surf * SA_surface/SA_rab11;
        parameters.kR2vegfA_R2on_surf = 0;
        parameters.kR2vegfA_R2on_rab4 = parameters.kR2vegfA_R2on_surf * SA_surface/SA_rab4;
        parameters.kR2vegfA_R2on_rab11 = parameters.kR2vegfA_R2on_surf * SA_surface/SA_rab11;
        parameters.kN1vegfA_N1on_surf = 0;
        parameters.kN1vegfA_N1on_rab4 = parameters.kN1vegfA_N1on_surf * SA_surface/SA_rab4;
        parameters.kN1vegfA_N1on_rab11 = parameters.kN1vegfA_N1on_surf * SA_surface/SA_rab11;
        parameters.kN1vegfA_N1off = parameters.kvegfA_N1off;
        parameters.kN1plgf2_N1on_surf = 0;
        parameters.kN1plgf2_N1on_rab4 = parameters.kN1plgf2_N1on_surf * SA_surface/SA_rab4;
        parameters.kN1plgf2_N1on_rab11 = parameters.kN1plgf2_N1on_surf * SA_surface/SA_rab11;
        parameters.kN1plgf2_N1off = parameters.kplgf_N1off;

        % Permit RVN coupling 
        parameters.kN1vegfA_R2on_surf = parameters.kN1vegfA_R2on_surf * 2; 
        parameters.kN1vegfA_R2on_rab4 = parameters.kN1vegfA_R2on_surf * SA_surface/SA_rab4; %parameters.kR2vegfA_R2on_rab4;
        parameters.kN1vegfA_R2on_rab11 = parameters.kN1vegfA_R2on_surf * SA_surface/SA_rab4; %parameters.kR2vegfA_R2on_rab11;
        parameters.kN1vegfA_R2off = 1e-3;

        parameters.kR2vegfA_N1on_surf = parameters.kR2vegfA_N1on_surf * 2; 
        parameters.kR2vegfA_N1on_rab4 = parameters.kR2vegfA_N1on_surf * SA_surface/SA_rab4;
        parameters.kR2vegfA_N1on_rab11 = parameters.kR2vegfA_N1on_surf * SA_surface/SA_rab11;
        parameters.kR2vegfA_N1off = 1e-3;

        % remove internal coupling
        parameters.kdeltaR1R1_V = 0;
        parameters.kdeltaR1R1_P = 0;
        parameters.kdeltaR2R2 = 0;
        parameters.kdeltaVEGFA_R = 0;
        parameters.kdeltaPLGF_R = 0;

end

%% RUN NO-LIGAND STEADY STATE

% STEADY STATE - ALL RECEPTORS PRESENT (unless excluded by tag above)
disp('running presimulation of receptors to steady state')

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

%% DEFINE PARAMETER SCENARIOS TO RUN

% LIGAND DOSE
v165a = 0; % [6843*1000*ngml./50];
v121a = 0; % [10750*1000*ngml./50]; 
plgf1 = 0; % [10140*1000*ngml./50]; 
plgf2 = 0; % [8702*1000*ngml./50];
switch ligandforthisrun 
    case "v165a"
        v165a = [6843*1000*ngml./50];
    case "v121a"
        v121a = [10750*1000*ngml./50]; 
    case "plgf1"
        plgf1 = [10140*1000*ngml./50]; 
    case "plgf2";
        plgf2 = [8702*1000*ngml./50];
end

% TRAFFICKING ADJUSTMENTS - SET UP ALL RUNS
traffnames = {'int','deg','rec4','rec4t11','rec11'};
ligandjnames = {'R1vegf','R1plgf','R2vegf','N1vegf','N1plgf'};
adjparamnames = {};
for i = 1:length(ligandjnames) % VR1, PR1, VR2, VN1, PN1
    for j = 1:length(traffnames) % int, deg, rec4, 4to11, rec11
        adjparamnames{length(traffnames)*(i-1)+j} = strcat('adjk',traffnames{j},ligandjnames{i});
    end
end

totalruns = length(traffnames)*length(ligandjnames)*length(paramadj);
for i = 1:(length(adjparamnames))
    scanparam.(adjparamnames{i})= 1.0; % initialize all with no change (1.0)
end
scanparams(1:totalruns) = scanparam; % set up a structure for each run

for i = 1:length(adjparamnames) % each of the 25 trafficking parameter adjustments
        for k = 1:length(paramadj) % values of adjustment
            scanparams(5*(i-1)+k).(adjparamnames{i}) = paramadj(k);
            if i == 11 % kintR2vegf
                scanparams(5*(i-1)+k).(adjparamnames{i}) = paramadjkintR2(k); % kintR2 is the only different set
            end
        end
end

% RUN THE CASES
disp('running trafficking cases')
for i = 1:totalruns
    runparams = adjustligandtraffparameters(parameters,scanparams(i));
    Routs(i) = RunLigandCase(model,speciesInit,runparams,v121a,v165a,plgf1,plgf2,timepts);
    disp(strcat('Trafficking case simulations =',num2str(i*100/totalruns),'% complete'))
    [R1rates(i),R2rates(i),N1rates(i)] = CalcRatesL(Routs(i),runparams); % perturbs?
    Routsnorm(i) = CalcNormOutputs(Routs(i));
end

outnames = fieldnames(Routs);
ratenames = fieldnames(R1rates);

mapr = [linspace(0.33,1,101)  linspace(.99,0,100)];
mapg = [linspace(0,1,101)   linspace(.99,.33,100)];
mapb = [linspace(0.33,1,101) linspace(0.99,0,100)];
map =[mapr' mapg' mapb'];


%% OUTPUT TIMECOURSES FOR THE METRICS

if longfigs == 1
    
disp('creating figures')

figoutnames = {'VEGFR1','VEGFR2','NRP1'};
paneloutnames.VEGFR1 = {'R1_surf','R1_int','R1_total'};
paneloutnames.VEGFR2 = {'R2_surf','R2_int','R2_total'};
paneloutnames.NRP1   = {'N1_surf','N1_int','N1_total'};

% Nine receptor outputs vs timepoints (for each parameter change) -
% absolute numbers
for i = 1:length(adjparamnames)
    for j = 1:length(figoutnames)
        f1 = genfigure(i, j, adjparamnames, figoutnames, paneloutnames, paramadj, timepts, Routs, "abs", scanparams, visibility);
    end
end

% Nine receptor outputs vs timepoints (for each parameter change) -
% numbers as a percentage of control
for i = 1:length(adjparamnames)
    for j = 1:length(figoutnames)
        f1 = genfigure(i, j, adjparamnames, figoutnames, paneloutnames, paramadj, timepts, Routsnorm, "norm", scanparams, visibility);
    end
end


end


%% NEXT - flux with x-axis being time, and three panels (R1, R2, N1) - one fig for each ligand/dose? (50 ng/ml)

%% ALSO - compare the fluxes to the changes in levels over time (never really showing the N1 yet?)

%% THEN - FUNCTIONALIZE these graphs so we can just call specific ones on loop



%% Function - Run VEGF x PlGF simulation array

% function [Routs] = RunVEGFPlGFarray(runModel,speciesInit,parameters,lig1,lig2,lig1conc,lig2conc,timepts)
function [ConcOut] = RunLigandCase(runModel,inits,parameters,v121a,v165a,plgf1,plgf2,Timepoints)

inits(1) = inits(1) + v165a;
inits(2) = inits(2) + v121a;
inits(3) = inits(3) + plgf1;
inits(4) = inits(4) + plgf2;

%% RUN MODEL
% [timepoints, species_out, observables_out] = Ligand_Traff_model_14Sep2023(timeLen, timestp_sim, parameters, inits);
[~, timepoints, species_out, observables_out] = eval(append(runModel,"(Timepoints', inits, parameters, 1)"));

ConcOut=CalcOutputsL(observables_out);

numpercell_to_pM_surf  = (1.0e4)/(6.022)/1.0e7; % => (#/cell) * (1e12 pmol/mol) * (1e15 fL/L) /(Na #/mol) / (fL/cell) 
numpercell_to_pM_int   = (1.0e4)/(6.022)/15; 
numpercell_to_pM_rab4  = (1.0e4)/(6.022)/11.25; 
numpercell_to_pM_rab11 = (1.0e4)/(6.022)/3.75; 

ConcOut.V165_surf  = ConcOut.V165_surf * numpercell_to_pM_surf;
ConcOut.V165_int   = ConcOut.V165_int * numpercell_to_pM_int;
ConcOut.V165_Rab4  = ConcOut.V165_Rab4 * numpercell_to_pM_rab4;
ConcOut.V165_Rab11 = ConcOut.V165_Rab11 * numpercell_to_pM_rab11;
ConcOut.V121_surf  = ConcOut.V121_surf * numpercell_to_pM_surf;
ConcOut.V121_int   = ConcOut.V121_int * numpercell_to_pM_int;
ConcOut.V121_Rab4  = ConcOut.V121_Rab4 * numpercell_to_pM_rab4;
ConcOut.V121_Rab11 = ConcOut.V121_Rab11 * numpercell_to_pM_rab11;
ConcOut.PGF1_surf  = ConcOut.PGF1_surf * numpercell_to_pM_surf;
ConcOut.PGF1_int   = ConcOut.PGF1_int * numpercell_to_pM_int;
ConcOut.PGF1_Rab4  = ConcOut.PGF1_Rab4 * numpercell_to_pM_rab4;
ConcOut.PGF1_Rab11 = ConcOut.PGF1_Rab11 * numpercell_to_pM_rab11;
ConcOut.PGF2_surf  = ConcOut.PGF2_surf * numpercell_to_pM_surf;
ConcOut.PGF2_int   = ConcOut.PGF2_int * numpercell_to_pM_int;
ConcOut.PGF2_Rab4  = ConcOut.PGF2_Rab4 * numpercell_to_pM_rab4;
ConcOut.PGF2_Rab11 = ConcOut.PGF2_Rab11 * numpercell_to_pM_rab11;

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

%% Function - generate figure

function f1 = genfigure(i, j, adjparamnames, figoutnames, paneloutnames, paramadj, timepts, Routs, absornorm, scanparams, visibility)

    figname = strcat("fig_",absornorm,"_",num2str(j),"_",num2str(i),"_",figoutnames{j},"_",adjparamnames{i});
    eval(strcat(figname," = figure('visible',visibility);"));
    calcname = strcat("outputs_",figname);

    locale = ["surf", "int", "total"];
    for n = 1:3
        axname = strcat("ax",num2str(n));
        eval(strcat(axname," = subplot(1,3,",num2str(n),");")); %     ax1 = subplot(1,3,1);
        a = [];
        for k = 1:length(paramadj)
            m = 5*(i-1)+k;
            outputa = figoutnames{j};
            outputb = string(paneloutnames.(outputa)(n));
            a = [a Routs(m).(outputb)(:)];
        end
        plot(timepts/60,a,'LineWidth',2.0); % (:,1,:) vegf and timepoints; % (1,:,:) plgf and timepoints
        hold on;
        gca.FontSize = 14;
        title(eval(axname), strcat(figoutnames{j},' Output ',outputb));
        xlabel(eval(axname),'Time (min)');
        ylabel(eval(axname), outputb);
        paramtemp = [];
        for k = 1:length(paramadj)
            paramtemp = [paramtemp scanparams(5*(i-1)+k).(adjparamnames{i})];
        end
        a = [0 paramtemp;timepts'/60 a];
        csvname = strcat("outputDataDr1/",calcname,"_",locale(n),".csv");
        csvwrite(csvname,a);
    end
  
	set(eval(figname), 'Position',[0 0 1500 600])
    set(eval(figname), 'PaperUnits', 'inches');
    set(eval(figname), 'PaperSize', [4 6]);
    exportgraphics(eval(figname),strcat("outputFigsDr1/",calcname,".png"),'Resolution',300);
    disp(strcat('Visualizations (abs) =',num2str((length(figoutnames)*(i-1)+j)*100/(length(adjparamnames)*length(figoutnames))),'% complete'))
    f1 = eval(figname); 
end
