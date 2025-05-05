close all;
clear all;

model = "All_Lig_Traff_model_Bionetgen_20231016"; %"TraffPhosph"
parameters = base_parameters(model);

% HIGHER VEGFR2 INTERNALIZATION
parameters.kVegfR2Rab5a = parameters.kR2Rab5a * 3;
parameters.kVegfR2N1Rab5a = parameters.kVegfR2Rab5a;

% ngml = [0 0.01 0.025 0.05 0.1 0.25 0.5 1 2.5 5 10 25 50 100 150];
ngml = 50; % standard ligand initial concentration
ligands = ["v121a","v165a","plgf1","plgf2"];
% ligandforthisrun = "v165a"; % alternatives "v121a"; "plgf1"; "plgf2";

paramadj = 5; % for adjusting trafficking parameters

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


% TRAFFICKING ADJUSTMENTS - SET UP ALL RUNS
traffnames = {'int','deg','rec4','rec4t11','rec11'};
ligandjnames = {'R1vegf','R1plgf','R2vegf','N1vegf','N1plgf'};
adjparamnames = {};
for i = 1:length(ligandjnames) % VR1, PR1, VR2, VN1, PN1
    for j = 1:length(traffnames) % int, deg, rec4, 4to11, rec11
        adjparamnames{length(traffnames)*(i-1)+j} = strcat('adjk',traffnames{j},ligandjnames{i});
    end
end

totalruns = length(ligands)*(length(traffnames)*length(ligandjnames)+1); % added 1 for control case
for i = 1:(length(adjparamnames))
    scanparam.(adjparamnames{i})= 1.0; % initialize all with no change (1.0)
end
scanparams(1:totalruns) = scanparam; % set up a structure for each run

for k = 1:length(ligands)
    % skip 1st one (that's the baseline/control)
    for i = 1:length(adjparamnames) % each of the 25 trafficking parameter adjustments
        scanparams((k-1)*(length(traffnames)*length(ligandjnames)+1)+i+1).(adjparamnames{i}) = paramadj;
    end
end

% for n = 1:length(ligands)
%     ligandforthisrun = ligands(n);
% 
% LIGAND DOSE
v165a(1:totalruns) = 0; % [6843*1000*ngml./50];
v121a = v165a; % [10750*1000*ngml./50]; 
plgf1 = v165a; % [10140*1000*ngml./50]; 
plgf2 = v165a; % [8702*1000*ngml./50];

offset = (length(traffnames)*length(ligandjnames)+1);
for i = 1:offset
    v121a(i) = [10750*1000*ngml./50]; 
    v165a(i+offset) = [6843*1000*ngml./50];
    plgf1(i+2*offset) = [10140*1000*ngml./50]; 
    plgf2(i+3*offset) = [8702*1000*ngml./50];
end

% RUN THE CASES
disp('running trafficking cases')
for i = 1:totalruns
    runparams = adjustligandtraffparameters(parameters,scanparams(i));
    Routs(i) = RunLigandCase(model,speciesInit,runparams,v121a(i),v165a(i),plgf1(i),plgf2(i),timepts);
    disp(strcat('Trafficking case simulations =',num2str(i*100/totalruns),'% complete'))
    [R1rates(i),R2rates(i),N1rates(i)] = CalcRatesL(Routs(i),runparams); % perturbs?
    Routsnorm(i) = CalcNormOutputs(Routs(i));    
end

for i = 1:offset
    Routsnormctrl(i) = CalcNormCtrlOutputs(Routsnorm(i),Routsnorm(1));
    Routsnormctrl(i+offset) = CalcNormCtrlOutputs(Routsnorm(i+offset),Routsnorm(1+offset));
    Routsnormctrl(i+2*offset) = CalcNormCtrlOutputs(Routsnorm(i+2*offset),Routsnorm(1+2*offset));
    Routsnormctrl(i+3*offset) = CalcNormCtrlOutputs(Routsnorm(i+3*offset),Routsnorm(1+3*offset));
    Routsctrl(i) = CalcNormCtrlOutputs(Routs(i),Routs(1));
    Routsctrl(i+offset) = CalcNormCtrlOutputs(Routs(i+offset),Routs(1+offset));
    Routsctrl(i+2*offset) = CalcNormCtrlOutputs(Routs(i+2*offset),Routs(1+2*offset));
    Routsctrl(i+3*offset) = CalcNormCtrlOutputs(Routs(i+3*offset),Routs(1+3*offset));
end

outnames = fieldnames(Routs);
ratenames = fieldnames(R1rates);

% two-color divergent map (down = purple; up = green)
mapr = [linspace(0.33,1,101)  linspace(.99,0,100)];
mapg = [linspace(0,1,101)   linspace(.99,.33,100)];
mapb = [linspace(0.33,1,101) linspace(0.99,0,100)];
map =[mapr' mapg' mapb'];

% one-color scale map (low = white; high = aqua blue)
mapr = linspace(1,66/255,201);
mapg = linspace(1,149/255,201); 
mapb = linspace(1,1,201); 
map1 =[mapr' mapg' mapb'];

ligandmap = [241 156 156; 191 107 104; 164 203 250; 18 51 98]/255;
barligandmap = [241 156 156; 191 107 104; 241 156 156; 191 107 104; 164 203 250; 18 51 98]/255;
barligandmap2 = [0 0 0; 241 156 156; 191 107 104; 164 203 250; 18 51 98]/255;

%% OUTPUT TIMECOURSES FOR THE METRICS

if longfigs == 1
    
disp('creating timecourse figures')

% Bar Charts -- active receptors

a1 = zeros(length(timepts),length(ligands)); 
a2 = a1; b1 = a1; b2 = a1; c1 = a1; c2 = a2;

for i = 1:length(ligands)
    % get the timecourses from the control case (no trafficking perturbation)
    a1(:,i) = Routs((i-1)*(length(traffnames)*length(ligandjnames)+1)+1).R1_surf_ACT;
    b1(:,i) = Routs((i-1)*(length(traffnames)*length(ligandjnames)+1)+1).R1_int_ACT;
    c1(:,i) = Routs((i-1)*(length(traffnames)*length(ligandjnames)+1)+1).R1_total_ACT;

    a2(:,i) = Routs((i-1)*(length(traffnames)*length(ligandjnames)+1)+1).R2_surf_ACT;
    b2(:,i) = Routs((i-1)*(length(traffnames)*length(ligandjnames)+1)+1).R2_int_ACT;
    c2(:,i) = Routs((i-1)*(length(traffnames)*length(ligandjnames)+1)+1).R2_total_ACT;
    
    aL(:,i) = Routs((i-1)*(length(traffnames)*length(ligandjnames)+1)+1).ligand_surf/1000; % divide by 1000 to get nM instead of pM;
    bL(:,i) = Routs((i-1)*(length(traffnames)*length(ligandjnames)+1)+1).ligand_int/1000; % divide by 1000 to get nM instead of pM;
%     cL(:,i) = Routs((i-1)*(length(traffnames)*length(ligandjnames)+1)+1).xyz;
end

bar1hrsurf = [a2(60+1,1:2) a1(60+1,1:4)];
bar4hrsurf = [a2(240+1,1:2) a1(240+1,1:4)];
bar1hrint = [b2(60+1,1:2) b1(60+1,1:4)];
bar4hrint = [b2(240+1,1:2) b1(240+1,1:4)];

barR11hrsurf = [a1(60+1,1:4)];
barR14hrsurf = [a1(240+1,1:4)];
barR11hrint = [b1(60+1,1:4)];
barR14hrint = [b1(240+1,1:4)];

barR21hrsurf = [a2(60+1,1:4)];
barR24hrsurf = [a2(240+1,1:4)];
barR21hrint = [b2(60+1,1:4)];
barR24hrint = [b2(240+1,1:4)];

bar1hrsurfL = [aL(60+1,1:4)];
bar4hrsurfL = [aL(240+1,1:4)];
bar1hrintL  = [bL(60+1,1:4)];
bar4hrintL  = [bL(240+1,1:4)];

figname = strcat("ActiveR_1hr4hr_bar");
eval(strcat(figname," = figure('visible',visibility);"));
calcname = strcat("outputs_",figname);

maxy = max(max([bar1hrsurf bar4hrsurf bar1hrint bar4hrint]));
ax1 = subplot(1,2,1);
bar1=bar([bar1hrsurf bar1hrint],'FaceColor','flat'); 
hold on;
bar1.CData = [barligandmap;barligandmap];
ylim(ax1,[0 maxy*1.03]);
ylabel(ax1,'Active Receptors (#/cell)');
colororder(ligandmap);

ax2 = subplot(1,2,2);
bar2=bar([bar4hrsurf bar4hrint],'FaceColor','flat'); 
hold on;
bar2.CData = [barligandmap;barligandmap];
ylim(ax2,[0 maxy*1.03]);
ylabel(ax2,'Active Receptors (#/cell)');
colororder(ligandmap);
csvname = strcat("outputDataDr1b/",calcname,".csv");
alldata = [bar1hrsurf bar1hrint bar4hrsurf bar4hrint];
csvwrite(csvname,alldata);

set(eval(figname), 'Position',[0 0 800 250])
set(eval(figname), 'PaperUnits', 'inches');
set(eval(figname), 'PaperSize', [4 6]);
exportgraphics(eval(figname),strcat("outputFigsDr1b/",calcname,".png"),'Resolution',300);

figname = strcat("ActiveR1_1hr4hr_bar");
eval(strcat(figname," = figure('visible',visibility);"));
calcname = strcat("outputs_",figname);

maxy = max(max([bar1hrsurf bar4hrsurf bar1hrint bar4hrint]));
ax1 = subplot(1,2,1);
bar1=bar([barR11hrsurf barR11hrint],'FaceColor','flat'); 
hold on;
bar1.CData = [ligandmap;ligandmap];
ylim(ax1,[0 maxy*1.03]);
ylabel(ax1,'Active VEGFR1 (#/cell)');
colororder(ligandmap);

ax2 = subplot(1,2,2);
bar2=bar([barR14hrsurf barR14hrint],'FaceColor','flat'); 
hold on;
bar2.CData = [ligandmap;ligandmap];
ylim(ax2,[0 maxy*1.03]);
ylabel(ax2,'Active VEGFR1 (#/cell)');
colororder(ligandmap);
csvname = strcat("outputDataDr1b/",calcname,".csv");
alldata = [barR11hrsurf barR11hrint barR14hrsurf barR14hrint];
csvwrite(csvname,alldata);

set(eval(figname), 'Position',[0 0 800 250])
set(eval(figname), 'PaperUnits', 'inches');
set(eval(figname), 'PaperSize', [4 6]);
exportgraphics(eval(figname),strcat("outputFigsDr1b/",calcname,".png"),'Resolution',300);

figname = strcat("ActiveR2_1hr4hr_bar");
eval(strcat(figname," = figure('visible',visibility);"));
calcname = strcat("outputs_",figname);

maxy = max(max([bar1hrsurf bar4hrsurf bar1hrint bar4hrint]));
ax1 = subplot(1,2,1);
bar1=bar([barR21hrsurf barR21hrint],'FaceColor','flat'); 
hold on;
bar1.CData = [ligandmap;ligandmap];
ylim(ax1,[0 maxy*1.03]);
ylabel(ax1,'Active VEGFR2 (#/cell)');
colororder(ligandmap);

ax2 = subplot(1,2,2);
bar2=bar([barR24hrsurf barR24hrint],'FaceColor','flat'); 
hold on;
bar2.CData = [ligandmap;ligandmap];
ylim(ax2,[0 maxy*1.03]);
ylabel(ax2,'Active VEGFR2 (#/cell)');
colororder(ligandmap);
csvname = strcat("outputDataDr1b/",calcname,".csv");
alldata = [barR21hrsurf barR21hrint barR24hrsurf barR24hrint];
csvwrite(csvname,alldata);

set(eval(figname), 'Position',[0 0 800 250])
set(eval(figname), 'PaperUnits', 'inches');
set(eval(figname), 'PaperSize', [4 6]);
exportgraphics(eval(figname),strcat("outputFigsDr1b/",calcname,".png"),'Resolution',300);


figname = strcat("LigandConc_1hr4hr_bar");
eval(strcat(figname," = figure('visible',visibility);"));
calcname = strcat("outputs_",figname);

maxy = max(max([bar1hrsurfL bar4hrsurfL bar1hrintL bar4hrintL]));
ax1 = subplot(1,2,1);
bar1=bar([bar1hrsurfL bar1hrintL],'FaceColor','flat'); 
hold on;
bar1.CData = [ligandmap;ligandmap];
ylim(ax1,[0 maxy*1.03]);
ylabel(ax1,'Ligand Concentration (nM)');
colororder(ligandmap);

ax2 = subplot(1,2,2);
bar2=bar([bar4hrsurfL bar4hrintL],'FaceColor','flat'); 
hold on;
bar2.CData = [ligandmap;ligandmap];
ylim(ax2,[0 maxy*1.03]);
ylabel(ax2,'Ligand Concentration (nM)');
colororder(ligandmap);
csvname = strcat("outputDataDr1b/",calcname,".csv");
alldata = [bar1hrsurfL bar1hrintL bar4hrsurfL bar4hrintL];
csvwrite(csvname,alldata);

set(eval(figname), 'Position',[0 0 800 250])
set(eval(figname), 'PaperUnits', 'inches');
set(eval(figname), 'PaperSize', [4 6]);
exportgraphics(eval(figname),strcat("outputFigsDr1b/",calcname,".png"),'Resolution',300);

% timecourse graphs - active receptors

figname = strcat("timecourses_R1active_fig");
eval(strcat(figname," = figure('visible',visibility);"));
calcname = strcat("outputs_",figname);

maxy = max(max([a1 b1 c1]));
ax1 = subplot(1,3,1);
plot(timepts/60,a1,'LineWidth',2.0);
hold on;
colororder(ligandmap);
legend('V121','V165','PlGF1','PlGF2','Location','northeast');
gca.FontSize = 14;
title(ax1,'Surface Active VEGFR1');
xlabel(ax1,'Time (min)');
ylabel(ax1,'Active Receptors (#/cell)');
xlim(ax1,[0 240]);
ylim(ax1,[0 maxy]);
csvname = strcat("outputDataDr1b/",calcname,"_a.csv");
csvwrite(csvname,a1);
  ax2 = subplot(1,3,2);
  plot(timepts/60,b1,'LineWidth',2.0); 
  hold on;
  colororder(ligandmap);
%   legend('V121','V165','PlGF1','PlGF2','Location','northwest');
  gca.FontSize = 14;
  title(ax2,'Internal Active VEGFR1');
  xlabel(ax2,'Time (min)');
  ylabel(ax2,'Active Receptors (#/cell)');
  xlim(ax2,[0 240]);
  ylim(ax2,[0 maxy]);
  csvname = strcat("outputDataDr1b/",calcname,"_b.csv");
  csvwrite(csvname,b1);
    ax3 = subplot(1,3,3);
    plot(timepts/60,c1,'LineWidth',2.0); 
    hold on;
    colororder(ligandmap);
%     legend('V121','V165','PlGF1','PlGF2','Location','northwest');
    gca.FontSize = 14;
    title(ax3,'Total Active VEGFR1');
    xlabel(ax3,'Time (min)');
    ylabel(ax3,'Active Receptors (#/cell)');
    xlim(ax3,[0 240]);
    ylim(ax3,[0 maxy]);
    csvname = strcat("outputDataDr1b/",calcname,"_c.csv");
    csvwrite(csvname,c1);
set(eval(figname), 'Position',[0 0 900 250])
set(eval(figname), 'PaperUnits', 'inches');
set(eval(figname), 'PaperSize', [4 6]);
exportgraphics(eval(figname),strcat("outputFigsDr1b/",calcname,".png"),'Resolution',300);

figname = strcat("timecourses_R2active_fig");
eval(strcat(figname," = figure('visible',visibility);"));
calcname = strcat("outputs_",figname);

maxy = max(max([a2 b2 c2]));
ax1 = subplot(1,3,1);
plot(timepts/60,a2,'LineWidth',2.0); 
hold on;
colororder(ligandmap);
legend('V121','V165','PlGF1','PlGF2','Location','northeast');
gca.FontSize = 14;
title(ax1,'Surface Active VEGFR2');
xlabel(ax1,'Time (min)');
ylabel(ax1,'Active Receptors (#/cell)');
xlim(ax1,[0 240]);
ylim(ax1,[0 maxy]);
csvname = strcat("outputDataDr1b/",calcname,"_a.csv");
csvwrite(csvname,a2);
  ax2 = subplot(1,3,2);
  plot(timepts/60,b2,'LineWidth',2.0); 
  hold on;
  colororder(ligandmap);
%   legend('V121','V165','PlGF1','PlGF2','Location','best');
  gca.FontSize = 14;
  title(ax2,'Internal Active VEGFR2');
  xlabel(ax2,'Time (min)');
  ylabel(ax2,'Active Receptors (#/cell)');
  xlim(ax2,[0 240]);
  ylim(ax2,[0 maxy]);
  csvname = strcat("outputDataDr1b/",calcname,"_b.csv");
  csvwrite(csvname,b2);
    ax3 = subplot(1,3,3);
    plot(timepts/60,c2,'LineWidth',2.0); 
    hold on;
    colororder(ligandmap);
%     legend('V121','V165','PlGF1','PlGF2','Location','best');
    gca.FontSize = 14;
    title(ax3,'Total Active VEGFR2');
    xlabel(ax3,'Time (min)');
    ylabel(ax3,'Active Receptors (#/cell)');
    xlim(ax3,[0 240]);
    ylim(ax3,[0 maxy]);
    csvname = strcat("outputDataDr1b/",calcname,"_c.csv");
    csvwrite(csvname,c2);
set(eval(figname), 'Position',[0 0 900 250])
set(eval(figname), 'PaperUnits', 'inches');
set(eval(figname), 'PaperSize', [4 6]);
exportgraphics(eval(figname),strcat("outputFigsDr1b/",calcname,".png"),'Resolution',300);

figname = strcat("timecourses_ligandconc_fig");
eval(strcat(figname," = figure('visible',visibility);"));
calcname = strcat("outputs_",figname);

maxy = max(max([aL bL]));
ax1 = subplot(1,2,1);
plot(timepts/60,aL,'LineWidth',2.0); 
hold on;
colororder(ligandmap);
legend('V121','V165','PlGF1','PlGF2','Location','southeast');
gca.FontSize = 14;
title(ax1,'Extracellular Ligand');
xlabel(ax1,'Time (min)');
ylabel(ax1,'Ligand Concentration (nM)'); 
xlim(ax1,[0 240]);
ylim(ax1,[0 maxy]);
csvname = strcat("outputDataDr1b/",calcname,"_a.csv");
csvwrite(csvname,a2);
  ax2 = subplot(1,2,2);
  plot(timepts/60,bL,'LineWidth',2.0); 
  hold on;
  colororder(ligandmap);
  legend('V121','V165','PlGF1','PlGF2','Location','northeast');
  gca.FontSize = 14;
  title(ax2,'Internal Ligand');
  xlabel(ax2,'Time (min)');
  ylabel(ax2,'Ligand Concentration (nM)');
  xlim(ax2,[0 240]);
  ylim(ax2,[0 maxy]);
  csvname = strcat("outputDataDr1b/",calcname,"_b.csv");
  csvwrite(csvname,b2);
set(eval(figname), 'Position',[0 0 600 250])
set(eval(figname), 'PaperUnits', 'inches');
set(eval(figname), 'PaperSize', [4 6]);
exportgraphics(eval(figname),strcat("outputFigsDr1b/",calcname,".png"),'Resolution',300);


% Bar Charts -- total receptors

a1 = zeros(length(timepts),length(ligands)); 
a2 = a1; b1 = a1; b2 = a1; c1 = a1; c2 = a2;

for i = 1:length(ligands)
    % get the timecourses from the control case (no trafficking perturbation)
    a1(:,i) = Routs((i-1)*(length(traffnames)*length(ligandjnames)+1)+1).R1_surf;
    b1(:,i) = Routs((i-1)*(length(traffnames)*length(ligandjnames)+1)+1).R1_int;
    c1(:,i) = Routs((i-1)*(length(traffnames)*length(ligandjnames)+1)+1).R1_total;

    a2(:,i) = Routs((i-1)*(length(traffnames)*length(ligandjnames)+1)+1).R2_surf;
    b2(:,i) = Routs((i-1)*(length(traffnames)*length(ligandjnames)+1)+1).R2_int;
    c2(:,i) = Routs((i-1)*(length(traffnames)*length(ligandjnames)+1)+1).R2_total;

    a3(:,i) = Routs((i-1)*(length(traffnames)*length(ligandjnames)+1)+1).N1_surf;
    b3(:,i) = Routs((i-1)*(length(traffnames)*length(ligandjnames)+1)+1).N1_int;
    c3(:,i) = Routs((i-1)*(length(traffnames)*length(ligandjnames)+1)+1).N1_total;
end

bar1hrsurf = [a2(60+1,1:2) a1(60+1,1:4)];
bar4hrsurf = [a2(240+1,1:2) a1(240+1,1:4)];
bar1hrint = [b2(60+1,1:2) b1(60+1,1:4)];
bar4hrint = [b2(240+1,1:2) b1(240+1,1:4)];

figname = strcat("TotalR_1hr4hr_bar");
eval(strcat(figname," = figure('visible',visibility);"));
calcname = strcat("outputs_",figname);

maxy = max(max([bar1hrsurf bar4hrsurf bar1hrint bar4hrint]));
ax1 = subplot(1,2,1);
bar1=bar([bar1hrsurf bar1hrint],'FaceColor','flat'); 
hold on;
bar1.CData = [barligandmap;barligandmap];
ylim(ax1,[0 maxy*1.03]);
ylabel(ax1,'Receptors (#/cell)');
colororder(ligandmap);

ax2 = subplot(1,2,2);
bar2=bar([bar4hrsurf bar4hrint],'FaceColor','flat'); 
hold on;
bar2.CData = [barligandmap;barligandmap];
ylim(ax2,[0 maxy*1.03]);
ylabel(ax2,'Receptors (#/cell)');
colororder(ligandmap);
csvname = strcat("outputDataDr1b/",calcname,".csv");
alldata = [bar1hrsurf bar1hrint bar4hrsurf bar4hrint];
csvwrite(csvname,alldata);

set(eval(figname), 'Position',[0 0 800 250])
set(eval(figname), 'PaperUnits', 'inches');
set(eval(figname), 'PaperSize', [4 6]);
exportgraphics(eval(figname),strcat("outputFigsDr1b/",calcname,".png"),'Resolution',300);

% timecourse graphs - total receptors

figname = strcat("timecourses_R1_fig");
eval(strcat(figname," = figure('visible',visibility);"));
calcname = strcat("outputs_",figname);

maxy = max(max([a1 b1 c1]));
ax1 = subplot(1,3,1);
plot(timepts/60,a1,'LineWidth',2.0);
hold on;
colororder(ligandmap);
legend('V121','V165','PlGF1','PlGF2','Location','northeast');
gca.FontSize = 14;
title(ax1,'Surface VEGFR1');
xlabel(ax1,'Time (min)');
ylabel(ax1,'Receptors (#/cell)');
xlim(ax1,[0 240]);
ylim(ax1,[0 maxy]);
csvname = strcat("outputDataDr1b/",calcname,"_a.csv");
csvwrite(csvname,a1);
  ax2 = subplot(1,3,2);
  plot(timepts/60,b1,'LineWidth',2.0); 
  hold on;
  colororder(ligandmap);
%   legend('V121','V165','PlGF1','PlGF2','Location','northwest');
  gca.FontSize = 14;
  title(ax2,'Internal VEGFR1');
  xlabel(ax2,'Time (min)');
  ylabel(ax2,'Receptors (#/cell)');
  xlim(ax2,[0 240]);
  ylim(ax2,[0 maxy]);
  csvname = strcat("outputDataDr1b/",calcname,"_b.csv");
  csvwrite(csvname,b1);
    ax3 = subplot(1,3,3);
    plot(timepts/60,c1,'LineWidth',2.0); 
    hold on;
    colororder(ligandmap);
%     legend('V121','V165','PlGF1','PlGF2','Location','northwest');
    gca.FontSize = 14;
    title(ax3,'Total VEGFR1');
    xlabel(ax3,'Time (min)');
    ylabel(ax3,'Receptors (#/cell)');
    xlim(ax3,[0 240]);
    ylim(ax3,[0 maxy]);
    csvname = strcat("outputDataDr1b/",calcname,"_c.csv");
    csvwrite(csvname,c1);
set(eval(figname), 'Position',[0 0 900 250])
set(eval(figname), 'PaperUnits', 'inches');
set(eval(figname), 'PaperSize', [4 6]);
exportgraphics(eval(figname),strcat("outputFigsDr1b/",calcname,".png"),'Resolution',300);

figname = strcat("timecourses_R2_fig");
eval(strcat(figname," = figure('visible',visibility);"));
calcname = strcat("outputs_",figname);

maxy = max(max([a2 b2 c2]));
ax1 = subplot(1,3,1);
plot(timepts/60,a2,'LineWidth',2.0); 
hold on;
colororder(ligandmap);
legend('V121','V165','PlGF1','PlGF2','Location','northeast');
gca.FontSize = 14;
title(ax1,'Surface VEGFR2');
xlabel(ax1,'Time (min)');
ylabel(ax1,'Receptors (#/cell)');
xlim(ax1,[0 240]);
ylim(ax1,[0 maxy]);
csvname = strcat("outputDataDr1b/",calcname,"_a.csv");
csvwrite(csvname,a2);
  ax2 = subplot(1,3,2);
  plot(timepts/60,b2,'LineWidth',2.0); 
  hold on;
  colororder(ligandmap);
%   legend('V121','V165','PlGF1','PlGF2','Location','best');
  gca.FontSize = 14;
  title(ax2,'Internal VEGFR2');
  xlabel(ax2,'Time (min)');
  ylabel(ax2,'Receptors (#/cell)');
  xlim(ax2,[0 240]);
  ylim(ax2,[0 maxy]);
  csvname = strcat("outputDataDr1b/",calcname,"_b.csv");
  csvwrite(csvname,b2);
    ax3 = subplot(1,3,3);
    plot(timepts/60,c2,'LineWidth',2.0); 
    hold on;
    colororder(ligandmap);
%     legend('V121','V165','PlGF1','PlGF2','Location','best');
    gca.FontSize = 14;
    title(ax3,'Total VEGFR2');
    xlabel(ax3,'Time (min)');
    ylabel(ax3,'Receptors (#/cell)');
    xlim(ax3,[0 240]);
    ylim(ax3,[0 maxy]);
    csvname = strcat("outputDataDr1b/",calcname,"_c.csv");
    csvwrite(csvname,c2);
set(eval(figname), 'Position',[0 0 900 250])
set(eval(figname), 'PaperUnits', 'inches');
set(eval(figname), 'PaperSize', [4 6]);
exportgraphics(eval(figname),strcat("outputFigsDr1b/",calcname,".png"),'Resolution',300);

figname = strcat("timecourses_N1_fig");
eval(strcat(figname," = figure('visible',visibility);"));
calcname = strcat("outputs_",figname);

maxy = max(max([a3 b3 c3]));
ax1 = subplot(1,3,1);
plot(timepts/60,a3,'LineWidth',2.0); 
hold on;
colororder(ligandmap);
legend('V121','V165','PlGF1','PlGF2','Location','northeast');
gca.FontSize = 14;
title(ax1,'Surface NRP1');
xlabel(ax1,'Time (min)');
ylabel(ax1,'Receptors (#/cell)');
xlim(ax1,[0 240]);
ylim(ax1,[0 maxy]);
csvname = strcat("outputDataDr1b/",calcname,"_a.csv");
csvwrite(csvname,a3);
  ax2 = subplot(1,3,2);
  plot(timepts/60,b3,'LineWidth',2.0); 
  hold on;
  colororder(ligandmap);
%   legend('V121','V165','PlGF1','PlGF2','Location','best');
  gca.FontSize = 14;
  title(ax2,'Internal NRP1');
  xlabel(ax2,'Time (min)');
  ylabel(ax2,'Receptors (#/cell)');
  xlim(ax2,[0 240]);
  ylim(ax2,[0 maxy]);
  csvname = strcat("outputDataDr1b/",calcname,"_b.csv");
  csvwrite(csvname,b3);
    ax3 = subplot(1,3,3);
    plot(timepts/60,c3,'LineWidth',2.0); 
    hold on;
    colororder(ligandmap);
%     legend('V121','V165','PlGF1','PlGF2','Location','best');
    gca.FontSize = 14;
    title(ax3,'Total NRP1');
    xlabel(ax3,'Time (min)');
    ylabel(ax3,'Receptors (#/cell)');
    xlim(ax3,[0 240]);
    ylim(ax3,[0 maxy]);
    csvname = strcat("outputDataDr1b/",calcname,"_c.csv");
    csvwrite(csvname,c3);
set(eval(figname), 'Position',[0 0 900 250])
set(eval(figname), 'PaperUnits', 'inches');
set(eval(figname), 'PaperSize', [4 6]);
exportgraphics(eval(figname),strcat("outputFigsDr1b/",calcname,".png"),'Resolution',300);



% Bar Charts -- surf-int receptor ratios

for i = 1:length(ligands)
    % get the timecourses from the control case (no trafficking perturbation)
    a1(:,i) = Routs((i-1)*(length(traffnames)*length(ligandjnames)+1)+1).R1_surf;
    b1(:,i) = Routs((i-1)*(length(traffnames)*length(ligandjnames)+1)+1).R2_surf;
    c1(:,i) = Routs((i-1)*(length(traffnames)*length(ligandjnames)+1)+1).N1_surf;

    a2(:,i) = Routs((i-1)*(length(traffnames)*length(ligandjnames)+1)+1).R1_int;
    b2(:,i) = Routs((i-1)*(length(traffnames)*length(ligandjnames)+1)+1).R2_int;
    c2(:,i) = Routs((i-1)*(length(traffnames)*length(ligandjnames)+1)+1).N1_int;
end

aratio = a1./a2; bratio = b1./b2; cratio = c1./c2;

barR1 = [aratio(1,1) aratio(60+1,1:4) aratio(1,1) aratio(240+1,1:4)]; % first one is control (no ligand)
barR2 = [bratio(1,1) bratio(60+1,1:4) bratio(1,1) bratio(240+1,1:4)];
barN1 = [cratio(1,1) cratio(60+1,1:4) cratio(1,1) cratio(240+1,1:4)];

barR1norm = barR1/barR1(1)*100;
barR2norm = barR2/barR2(1)*100;
barN1norm = barN1/barN1(1)*100;

figname = strcat("SurfIntRatio_bar");
eval(strcat(figname," = figure('visible',visibility);"));
calcname = strcat("outputs_",figname);

maxy = max(max([barR1 barR2 barN1]));
ax1 = subplot(1,3,1);
bar1=bar(barR1,'FaceColor','flat'); 
hold on;
bar1.CData = [barligandmap2;barligandmap2];
ylim(ax1,[0 maxy*1.03]);
ylabel(ax1,'Surface:Internal Ratio');
colororder(ligandmap);

ax2 = subplot(1,3,2);
bar2=bar(barR2,'FaceColor','flat'); 
hold on;
bar2.CData = [barligandmap2;barligandmap2];
ylim(ax2,[0 maxy*1.03]);
ylabel(ax2,'Surface:Internal Ratio');
colororder(ligandmap);

ax3 = subplot(1,3,3);
bar3=bar(barN1,'FaceColor','flat'); 
hold on;
bar3.CData = [barligandmap2;barligandmap2];
ylim(ax3,[0 maxy*1.03]);
ylabel(ax3,'Surface:Internal Ratio');
colororder(ligandmap);

csvname = strcat("outputDataDr1b/",calcname,".csv");
alldata = [barR1 barR2 barN1];
csvwrite(csvname,alldata);

set(eval(figname), 'Position',[0 0 800 250])
set(eval(figname), 'PaperUnits', 'inches');
set(eval(figname), 'PaperSize', [4 6]);
exportgraphics(eval(figname),strcat("outputFigsDr1b/",calcname,".png"),'Resolution',300);


figname = strcat("SurfIntRatioNorm_bar");
eval(strcat(figname," = figure('visible',visibility);"));
calcname = strcat("outputs_",figname);

maxy = max(max([barR1norm barR2norm barN1norm]));
ax1 = subplot(1,3,1);
bar1=bar(barR1norm,'FaceColor','flat'); 
hold on;
bar1.CData = [barligandmap2;barligandmap2];
ylim(ax1,[0 maxy*1.03]);
ylabel(ax1,'Surface:Internal Ratio (normalized)');
colororder(ligandmap);

ax2 = subplot(1,3,2);
bar2=bar(barR2norm,'FaceColor','flat'); 
hold on;
bar2.CData = [barligandmap2;barligandmap2];
ylim(ax2,[0 maxy*1.03]);
ylabel(ax2,'Surface:Internal Ratio (normalized)');
colororder(ligandmap);

ax3 = subplot(1,3,3);
bar3=bar(barN1norm,'FaceColor','flat'); 
hold on;
bar3.CData = [barligandmap2;barligandmap2];
ylim(ax3,[0 maxy*1.03]);
ylabel(ax3,'Surface:Internal Ratio (normalized)');
colororder(ligandmap);

csvname = strcat("outputDataDr1b/",calcname,".csv");
alldata = [barR1norm barR2norm barN1norm];
csvwrite(csvname,alldata);

set(eval(figname), 'Position',[0 0 800 250])
set(eval(figname), 'PaperUnits', 'inches');
set(eval(figname), 'PaperSize', [4 6]);
exportgraphics(eval(figname),strcat("outputFigsDr1b/",calcname,".png"),'Resolution',300);


end


%% Heatmaps

outputlist = ["R1_surf","R1_int","R1_total","R2_surf","R2_int","R2_total","N1_surf","N1_int","N1_total"];
times = [15, 60, 240];

% Nine receptor outputs vs timepoints (for each parameter change) -
% both relative numbers (compared to zero timepoint) 
% and relative numbers (compared to control)
for j = 1:length(outputlist)
    for k = 1:length(times)
        outputsxtimes((j-1)*length(times)+k) = strcat(outputlist(j),num2str(times(k)));
    end
end
for i = 1:length(ligands)
    heatmapmatrix1 = zeros((length(traffnames)*length(ligandjnames)+1),(length(outputlist)*length(times)));
    heatmapmatrix2 = heatmapmatrix1;
    for n = 1:(length(traffnames)*length(ligandjnames)+1)
        for j = 1:length(outputlist)
            % gather the norm and ctrl outputs at 15, 60, 240 mins
            for k = 1:length(times)
                heatmapmatrix1(n,(j-1)*length(times)+k) = Routsnorm((i-1)*(length(traffnames)*length(ligandjnames)+1)+n).(outputlist(j))(times(k)+1);
                heatmapmatrix2(n,(j-1)*length(times)+k) = Routsnormctrl((i-1)*(length(traffnames)*length(ligandjnames)+1)+n).(outputlist(j))(times(k)+1);
            end
        end
    end
% relative numbers (compared to zero timepoint) 
    f1 = figure;
    h1 = heatmap(outputsxtimes,['ctrl',adjparamnames],heatmapmatrix1); 
    h1.Colormap = map;    
	h1.ColorLimits = [0,200];
    h1.CellLabelColor = 'None';
    h1.MissingDataColor = [0.9 0.9 0.9];
    h1.MissingDataLabel = "n/a";
    h1.FontSize = 14;
%     h1.Title = 'Outputs - R1/R2/N1';
    h1.XLabel = 'outputs';
    h1.YLabel = 'parameters';
    set(f1, 'Position',[0 0 600 600])
    set(f1, 'PaperUnits', 'inches');
    set(f1, 'PaperSize', [4 4]);
    exportgraphics(f1,strcat("outputFigsDr1b/","RecLevels_",ligands(i),"_norm.png"),'Resolution',300);
%     a = [0 outputsxtimes; ['ctrl' convertCharsToStrings(adjparamnames)]' heatmapmatrix1];
    csvname = strcat("outputDataDr1b/","RecLevels_",ligands(i),"_norm.csv");
    csvwrite(csvname,heatmapmatrix1);
% relative numbers (compared to control)
    f2 = figure;
    h1 = heatmap(outputsxtimes,['ctrl',adjparamnames],heatmapmatrix2); 
    h1.Colormap = map;    
	h1.ColorLimits = [0,200];
    h1.CellLabelColor = 'None';
    h1.MissingDataColor = [0.9 0.9 0.9];
    h1.MissingDataLabel = "n/a";
    h1.FontSize = 14;
%     h1.Title = 'Outputs - R1/R2/N1';
    h1.XLabel = 'outputs';
    h1.YLabel = 'parameters';
    set(f2, 'Position',[0 0 600 600])
    set(f2, 'PaperUnits', 'inches');
    set(f2, 'PaperSize', [4 4]);
    exportgraphics(f2,strcat("outputFigsDr1b/","RecLevels_",ligands(i),"_normctrl.png"),'Resolution',300);
%     a = [0 outputsxtimes; ['ctrl' convertCharsToStrings(adjparamnames)]' heatmapmatrix2];
    csvname = strcat("outputDataDr1b/","RecLevels_",ligands(i),"_normctrl.csv");
    csvwrite(csvname,heatmapmatrix2);
end

outputlist = ["R1_surf_ACT", "R1_int_ACT", "R1_total_ACT"];

% Three activated receptor (R1_act) outputs vs timepoints (for each parameter change) -
% both absolute and relative numbers (compared to control)
for j = 1:length(outputlist)
    for k = 1:length(times)
        outputsxtimes2((j-1)*length(times)+k) = strcat(outputlist(j),num2str(times(k)));
    end
end
for i = 1:length(ligands)
    heatmapmatrix1 = zeros((length(traffnames)*length(ligandjnames)+1),(length(outputlist)*length(times)));
    heatmapmatrix2 = heatmapmatrix1;
    for n = 1:(length(traffnames)*length(ligandjnames)+1)
        for j = 1:length(outputlist)
            % gather the norm and ctrl outputs at 15, 60, 240 mins
            for k = 1:length(times)
                heatmapmatrix1(n,(j-1)*length(times)+k) = Routs((i-1)*(length(traffnames)*length(ligandjnames)+1)+n).(outputlist(j))(times(k)+1);
                heatmapmatrix2(n,(j-1)*length(times)+k) = Routsctrl((i-1)*(length(traffnames)*length(ligandjnames)+1)+n).(outputlist(j))(times(k)+1);
            end
        end
    end
% absolute numbers
    f1 = figure;
    h1 = heatmap(outputsxtimes2,['ctrl',adjparamnames],heatmapmatrix1); 
    h1.Colormap = map1;    
%     h1.ColorScaling = 'log';
%     lim1 = 0.1;
%     lim2 = max(1,ceil(log(max(max(heatmapmatrix1)))));
% %     [lim1,lim2]
%  	h1.ColorLimits = [0,lim2];
    lim2 = max(1,max(max(heatmapmatrix1)));
  	h1.ColorLimits = [0,lim2];
    h1.CellLabelColor = 'None';
    h1.MissingDataColor = [0.9 0.9 0.9];
    h1.MissingDataLabel = "n/a";
    h1.FontSize = 14;
%     h1.Title = 'Outputs - R1/R2/N1';
    h1.XLabel = 'outputs';
    h1.YLabel = 'parameters';
    set(f1, 'Position',[0 0 400 600])
    set(f1, 'PaperUnits', 'inches');
    set(f1, 'PaperSize', [4 4]);
    exportgraphics(f1,strcat("outputFigsDr1b/","ActiveRec_R1_",ligands(i),"_norm.png"),'Resolution',300);
%     b = [0 outputsxtimes2; ['ctrl' convertCharsToStrings(adjparamnames)]' heatmapmatrix1];
    csvname = strcat("outputDataDr1b/","ActiveRec_R1_",ligands(i),"_norm.csv");
    csvwrite(csvname,heatmapmatrix1);
% relative numbers (compared to control)
    f2 = figure;
    h1 = heatmap(outputsxtimes2,['ctrl',adjparamnames],heatmapmatrix2); 
    h1.Colormap = map;    
	h1.ColorLimits = [0,200];
    h1.CellLabelColor = 'None';
    h1.MissingDataColor = [0.9 0.9 0.9];
    h1.MissingDataLabel = "n/a";
    h1.FontSize = 14;
%     h1.Title = 'Outputs - R1/R2/N1';
    h1.XLabel = 'outputs';
    h1.YLabel = 'parameters';
    set(f2, 'Position',[0 0 400 600])
    set(f2, 'PaperUnits', 'inches');
    set(f2, 'PaperSize', [4 4]);
    exportgraphics(f2,strcat("outputFigsDr1b/","ActiveRec_R1_",ligands(i),"_normctrl.png"),'Resolution',300);
%     b = [0 outputsxtimes2; ['ctrl' convertCharsToStrings(adjparamnames)]' heatmapmatrix2];
    csvname = strcat("outputDataDr1b/","ActiveRec_R1_",ligands(i),"_normctrl.csv");
    csvwrite(csvname,heatmapmatrix2);
end

outputlist = ["R2_surf_ACT", "R2_int_ACT", "R2_total_ACT"];

% Three activated receptor (R2_act) outputs vs timepoints (for each parameter change) -
% both absolute and relative numbers
for j = 1:length(outputlist)
    for k = 1:length(times)
        outputsxtimes3((j-1)*length(times)+k) = strcat(outputlist(j),num2str(times(k)));
    end
end
for i = 1:length(ligands)
    heatmapmatrix1 = zeros((length(traffnames)*length(ligandjnames)+1),(length(outputlist)*length(times)));
    heatmapmatrix2 = heatmapmatrix1;
    for n = 1:(length(traffnames)*length(ligandjnames)+1)
        for j = 1:length(outputlist)
            % gather the norm and ctrl outputs at 15, 60, 240 mins
            for k = 1:length(times)
                heatmapmatrix1(n,(j-1)*length(times)+k) = Routs((i-1)*(length(traffnames)*length(ligandjnames)+1)+n).(outputlist(j))(times(k)+1);
                heatmapmatrix2(n,(j-1)*length(times)+k) = Routsctrl((i-1)*(length(traffnames)*length(ligandjnames)+1)+n).(outputlist(j))(times(k)+1);
            end
        end
    end
% absolute numbers
    f1 = figure;
    h1 = heatmap(outputsxtimes3,['ctrl',adjparamnames],heatmapmatrix1); 
    h1.Colormap = map1;    
%     h1.ColorScaling = 'log';
%     lim1 = 0.1;
%     lim2 = max(1,ceil(log(max(max(heatmapmatrix1)))));
% %     [lim1,lim2]
%  	h1.ColorLimits = [0,lim2];
    lim2 = max(1,max(max(heatmapmatrix1)));
  	h1.ColorLimits = [0,lim2];
    h1.CellLabelColor = 'None';
    h1.MissingDataColor = [0.9 0.9 0.9];
    h1.MissingDataLabel = "n/a";
    h1.FontSize = 14;
%     h1.Title = 'Outputs - R1/R2/N1';
    h1.XLabel = 'outputs';
    h1.YLabel = 'parameters';
    set(f1, 'Position',[0 0 400 600])
    set(f1, 'PaperUnits', 'inches');
    set(f1, 'PaperSize', [4 4]);
    exportgraphics(f1,strcat("outputFigsDr1b/","ActiveRec_R2_",ligands(i),"_norm.png"),'Resolution',300);
%     c = [0 outputsxtimes3; ['ctrl' convertCharsToStrings(adjparamnames)]' heatmapmatrix1];
    csvname = strcat("outputDataDr1b/","ActiveRec_R2_",ligands(i),"_norm.csv");
    csvwrite(csvname,heatmapmatrix1);
% relative numbers (compared to control)
    f2 = figure;
    h1 = heatmap(outputsxtimes3,['ctrl',adjparamnames],heatmapmatrix2); 
    h1.Colormap = map;    
	h1.ColorLimits = [0,200];
    h1.CellLabelColor = 'None';
    h1.MissingDataColor = [0.9 0.9 0.9];
    h1.MissingDataLabel = "n/a";
    h1.FontSize = 14;
%     h1.Title = 'Outputs - R1/R2/N1';
    h1.XLabel = 'outputs';
    h1.YLabel = 'parameters';
    set(f2, 'Position',[0 0 400 600])
    set(f2, 'PaperUnits', 'inches');
    set(f2, 'PaperSize', [4 4]);
    exportgraphics(f2,strcat("outputFigsDr1b/","ActiveRec_R2_",ligands(i),"_normctrl.png"),'Resolution',300);
%     c = [0 outputsxtimes3; ['ctrl' convertCharsToStrings(adjparamnames)]' heatmapmatrix2];
    csvname = strcat("outputDataDr1b/","ActiveRec_R2_",ligands(i),"_normctrl.csv");
    csvwrite(csvname,heatmapmatrix2);
end

% output the graphs!


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

ConcOut.ligand_surf = ConcOut.ligand_surf * numpercell_to_pM_surf;
ConcOut.ligand_int  = ConcOut.ligand_int  * numpercell_to_pM_int;

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

