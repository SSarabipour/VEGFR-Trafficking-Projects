close all;
clear all;

model = "All_Lig_Traff_model_Bionetgen_20231016"; %"TraffPhosph"
parameters = base_parameters(model);

% HIGHER VEGFR2 INTERNALIZATION
parameters.kVegfR2Rab5a = parameters.kR2Rab5a * 3;
parameters.kVegfR2N1Rab5a = parameters.kVegfR2Rab5a;

% ngml = [0 0.01 0.025 0.05 0.1 0.25 0.5 1 2.5 5 10 25 50 100 150];
ngml = [0 logspace(-2,2.3,100)];
minngml = min(ngml(2:end));
maxngml = max(ngml(2:end));
timepts = 0:60:(240*60); %[0 15 30 60 60 240];

ReFitProductionRates = 0 ; % 1 for yes, 0 for no
R1surf_target = 1800;
R2surf_target = 4900; 
N1surf_target = 68000;

visibility = 'off';
longfigs = 0; 

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

ScatchardMode = "off"; % when "on", assumes no production, no internalization, only binding (4dC)

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

% STEADY STATE - ALL RECEPTORS PRESENT
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

%% DEFINE LIGAND DOSE SCENARIOS TO RUN

% initial size
v165a = zeros(4*length(ngml),1);
v121a = v165a;
plgf1 = v165a;
plgf2 = v165a;

% dose range for each ligand, one at a time
span = 1:length(ngml);
v121a(span)                 = [10750*1000*ngml./50]; 
v165a(span+length(ngml))    = [6843*1000*ngml./50]; 
plgf1(span+2*length(ngml))  = [10140*1000*ngml./50]; 
plgf2(span+3*length(ngml))  = [8702*1000*ngml./50];

disp('running ligand cases')

for i = 1:length(v121a)
    Routs(i) = RunLigandCase(model,speciesInit,parameters,v121a(i),v165a(i),plgf1(i),plgf2(i),timepts);
    disp(strcat('Ligand cases =',num2str(i*100/(length(v121a))),'% complete'))
    [R1rates(i),R2rates(i),N1rates(i)] = CalcRatesL(Routs(i),parameters); % perturbs?
    Routsnorm(i) = CalcNormOutputs(Routs(i));
    R1ratesnorm(i) = CalcNormOutputs(R1rates(i));
    R2ratesnorm(i) = CalcNormOutputs(R2rates(i));
    N1ratesnorm(i) = CalcNormOutputs(N1rates(i));
end

outnames = fieldnames(Routs);
ratenames = fieldnames(R1rates);

mapr = [linspace(0.33,1,101)  linspace(.99,0,100)];
mapg = [linspace(0,1,101)   linspace(.99,.33,100)];
mapb = [linspace(0.33,1,101) linspace(0.99,0,100)];
map =[mapr' mapg' mapb'];

% blue,red,green maps for R1, R2, N1 lines
map1 = [linspace(.75,0,length(ngml))];
map2 = [linspace(1,1,length(ngml))];
mapblue =[map1' map1' map2'];
mapred =[map2' map1' map1'];
mapgreen =[map1' map2' map1'];
mapall = [mapblue;mapred;mapgreen];

ligandmap = [241 156 156; 191 107 104; 164 203 250; 18 51 98]/255;


%% OUTPUT TIMECOURSES FOR THE METRICS

if longfigs == 1
    
disp('creating figures')

% Output vs timepoints (for each ligand) (for each output)
for i = 1:length(outnames)
    figname = strcat("vafig",num2str(i));
    eval(strcat(figname," = figure('visible',visibility);"));
    calcname = strcat("ROUToverTime_fig_",num2str(i,'%02i'),"_",outnames{i});
    ax1 = subplot(4,1,1);
    a = [];
    for j = 1:length(ngml)
        a = [a Routs(j).(outnames{i})(:)];
    end
    plot(timepts/60,a); % (:,1,:) vegf and timepoints; % (1,:,:) plgf and timepoints
    hold on;
    gca.FontSize = 14;
    title(ax1, strcat('V121 - Output ',outnames{i}));
    xlabel(ax1,'Time (min)');
    ylabel(ax1, outnames{i});
	a = [0 ngml;timepts'/60 a];
    csvname = strcat("outputDataDr2/",calcname,"a.csv");
    csvwrite(csvname,a);
      ax2 = subplot(4,1,2);
      a = [];
      for j = ((1:length(ngml))+length(ngml))
        a = [a Routs(j).(outnames{i})(:)];
      end
      plot(timepts/60,a); % (:,1,:) vegf and timepoints; % (1,:,:) plgf and timepoints
      hold on;
      gca.FontSize = 14;
      title(ax2, strcat('V165 - Output ',outnames{i}));
      xlabel(ax2,'Time (min)');
      ylabel(ax2, outnames{i});
      a = [0 ngml;timepts'/60 a];
      csvname = strcat("outputDataDr2/",calcname,"b.csv");
      csvwrite(csvname,a);
        ax3 = subplot(4,1,3);
        a = [];
        for j = ((1:length(ngml))+2*length(ngml))
          a = [a Routs(j).(outnames{i})(:)];
        end
        plot(timepts/60,a); % (:,1,:) vegf and timepoints; % (1,:,:) plgf and timepoints
        hold on;
        gca.FontSize = 14;
        title(ax3, strcat('PlGF1 - Output ',outnames{i}));
        xlabel(ax3,'Time (min)');
        ylabel(ax3, outnames{i});
        a = [0 ngml;timepts'/60 a];
        csvname = strcat("outputDataDr2/",calcname,"c.csv");
        csvwrite(csvname,a);
          ax4 = subplot(4,1,4);
          a = [];
          for j = ((1:length(ngml))+3*length(ngml))
            a = [a Routs(j).(outnames{i})(:)];
          end
          plot(timepts/60,a); % (:,1,:) vegf and timepoints; % (1,:,:) plgf and timepoints
          hold on;
          gca.FontSize = 14;
          title(ax4, strcat('PlGF2 - Output ',outnames{i}));
          xlabel(ax4,'Time (min)');
          ylabel(ax4, outnames{i});
          a = [0 ngml;timepts'/60 a];
          csvname = strcat("outputDataDr2/",calcname,"d.csv");
          csvwrite(csvname,a);
    set(eval(figname), 'Position',[0 0 1000 1500])
    set(eval(figname), 'PaperUnits', 'inches');
    set(eval(figname), 'PaperSize', [4 6]);
    exportgraphics(eval(figname),strcat("outputFigsDr2/",calcname,".png"),'Resolution',300);
    disp(strcat('Visualizations =',num2str(i*100/(length(outnames))),'% complete'))
end

% Output vs timepoints (for each ligand) (for each output)
for i = 1:length(ratenames)
    figname = strcat("ratefig",num2str(i));
    eval(strcat(figname," = figure('visible',visibility);"));
    calcname = strcat("RatesoverTime_fig_",num2str(i,'%02i'),"_",ratenames{i});
    ax1 = subplot(4,1,1);
    a = []; b = []; c = [];
    for j = 1:length(ngml)
        a = [a R1rates(j).(ratenames{i})(:)];
        b = [b R2rates(j).(ratenames{i})(:)];
        c = [c N1rates(j).(ratenames{i})(:)];
    end
    plot(timepts/60,a); % (:,1,:) vegf and timepoints; % (1,:,:) plgf and timepoints
    hold on;
    plot(timepts/60,b); % (:,1,:) vegf and timepoints; % (1,:,:) plgf and timepoints
    plot(timepts/60,c); % (:,1,:) vegf and timepoints; % (1,:,:) plgf and timepoints
    colororder(mapall);
    gca.FontSize = 14;
    title(ax1, strcat('V121 - Output ',ratenames{i}));
    xlabel(ax1,'Time (min)');
    ylabel(ax1, ratenames{i});
	a = [0 ngml;timepts'/60 a];
	b = [0 ngml;timepts'/60 b];
	c = [0 ngml;timepts'/60 c];
    csvname = strcat("outputDataDr2/",calcname,"_R1_121.csv");
    csvwrite(csvname,a);
    csvname = strcat("outputDataDr2/",calcname,"_R2_121.csv");
    csvwrite(csvname,b);
    csvname = strcat("outputDataDr2/",calcname,"_N1_121.csv");
    csvwrite(csvname,c);
      ax2 = subplot(4,1,2);
      a = []; b = []; c = [];
      for j = ((1:length(ngml))+length(ngml))
        a = [a R1rates(j).(ratenames{i})(:)];
        b = [b R2rates(j).(ratenames{i})(:)];
        c = [c N1rates(j).(ratenames{i})(:)];
      end
      plot(timepts/60,a); % (:,1,:) vegf and timepoints; % (1,:,:) plgf and timepoints
      hold on;
      plot(timepts/60,b); % (:,1,:) vegf and timepoints; % (1,:,:) plgf and timepoints
      plot(timepts/60,c); % (:,1,:) vegf and timepoints; % (1,:,:) plgf and timepoints
      colororder(mapall);
      gca.FontSize = 14;
      title(ax2, strcat('V165 - Output ',ratenames{i}));
      xlabel(ax2,'Time (min)');
      ylabel(ax2, ratenames{i});
	a = [0 ngml;timepts'/60 a];
	b = [0 ngml;timepts'/60 b];
	c = [0 ngml;timepts'/60 c];
    csvname = strcat("outputDataDr2/",calcname,"_R1_165.csv");
    csvwrite(csvname,a);
    csvname = strcat("outputDataDr2/",calcname,"_R2_165.csv");
    csvwrite(csvname,b);
    csvname = strcat("outputDataDr2/",calcname,"_N1_165.csv");
    csvwrite(csvname,c);
        ax3 = subplot(4,1,3);
        a = []; b = []; c = [];
        for j = ((1:length(ngml))+2*length(ngml))
            a = [a R1rates(j).(ratenames{i})(:)];
            b = [b R2rates(j).(ratenames{i})(:)];
            c = [c N1rates(j).(ratenames{i})(:)];
        end
        plot(timepts/60,a); % (:,1,:) vegf and timepoints; % (1,:,:) plgf and timepoints
        hold on;
        plot(timepts/60,b); % (:,1,:) vegf and timepoints; % (1,:,:) plgf and timepoints
        plot(timepts/60,c); % (:,1,:) vegf and timepoints; % (1,:,:) plgf and timepoints
        colororder(mapall);
        gca.FontSize = 14;
        title(ax3, strcat('PlGF1 - Output ',ratenames{i}));
        xlabel(ax3,'Time (min)');
        ylabel(ax3, ratenames{i});
        a = [0 ngml;timepts'/60 a];
        b = [0 ngml;timepts'/60 b];
        c = [0 ngml;timepts'/60 c];
        csvname = strcat("outputDataDr2/",calcname,"_R1_P1.csv");
        csvwrite(csvname,a);
        csvname = strcat("outputDataDr2/",calcname,"_R2_P1.csv");
        csvwrite(csvname,b);
        csvname = strcat("outputDataDr2/",calcname,"_N1_P1.csv");
        csvwrite(csvname,c);
          ax4 = subplot(4,1,4);
          a = []; b = []; c = [];
          for j = ((1:length(ngml))+3*length(ngml))
            a = [a R1rates(j).(ratenames{i})(:)];
            b = [b R2rates(j).(ratenames{i})(:)];
            c = [c N1rates(j).(ratenames{i})(:)];
          end
          plot(timepts/60,a); % (:,1,:) vegf and timepoints; % (1,:,:) plgf and timepoints
          hold on;
          plot(timepts/60,b); % (:,1,:) vegf and timepoints; % (1,:,:) plgf and timepoints
          plot(timepts/60,c); % (:,1,:) vegf and timepoints; % (1,:,:) plgf and timepoints
          colororder(mapall);
          gca.FontSize = 14;
          title(ax4, strcat('PlGF2 - Output ',ratenames{i}));
          xlabel(ax4,'Time (min)');
          ylabel(ax4, ratenames{i});
        a = [0 ngml;timepts'/60 a];
        b = [0 ngml;timepts'/60 b];
        c = [0 ngml;timepts'/60 c];
        csvname = strcat("outputDataDr2/",calcname,"_R1_P2.csv");
        csvwrite(csvname,a);
        csvname = strcat("outputDataDr2/",calcname,"_R2_P2.csv");
        csvwrite(csvname,b);
        csvname = strcat("outputDataDr2/",calcname,"_N1_P2.csv");
        csvwrite(csvname,c);
    set(eval(figname), 'Position',[0 0 1000 1500])
    set(eval(figname), 'PaperUnits', 'inches');
    set(eval(figname), 'PaperSize', [4 6]);
    exportgraphics(eval(figname),strcat("outputFigsDr2/",calcname,".png"),'Resolution',300);
    disp(strcat('Visualizations =',num2str(i*100/(length(ratenames))),'% complete'))
end

end


%% Dimerization curves - total and surface

% R1, R2, N1 dimerization

RsDimznFig1 = VizRsdimerization(Routs,"R1_dimericfraction_surf","R2_dimericfraction_surf","N1_dimericfraction_surf","Dimzn_Rs_surf",ngml,visibility,"northwest",ligandmap);
RsDimznFig2 = VizRsdimerization(Routs,"R1_dimericfraction","R2_dimericfraction","N1_dimericfraction","Dimzn_Rs_tot",ngml,visibility,"northwest",ligandmap);

% R1-N1 Dimerization curves

R1N1DimznFig1 = VizR1N1dimerization(Routs,"R1_fractioninR1N1_surf","N1_fractioninR1N1_surf","Dimzn_R1N1_surf",ngml,visibility,"best",ligandmap);
R1N1DimznFig1 = VizR1N1dimerization(Routs,"R1_fractioninR1N1","N1_fractioninR1N1","Dimzn_R1N1_tot",ngml,visibility,"best",ligandmap);

%% Active receptor curves - by ligand as a check

activeFig1 = VizActiveV121V165P1P2(Routs,"R1R1_121a_total_ACT","R1R1_165a_total_ACT","R1R1_plgf1_total_ACT","R1R1_plgf2_total_ACT","Active R1 complexes","Active_R1_bylig",ngml,visibility,ligandmap);
activeFig2 = VizActiveV121V165(Routs,"R2R2_121a_total_ACT","R2R2_165a_total_ACT","Active R2 complexes","Active_R2_bylig",ngml,visibility,ligandmap);

%% R1, R2 active (total) 

activeFig3 = VizActiveV121V165P1P2(Routs,"R1_total_ACT","R1_total_ACT","R1_total_ACT","R1_total_ACT","Active R1 complexes","Active_R1_tot",ngml,visibility,ligandmap);
activeFig4 = VizActiveV121V165(Routs,"R2_total_ACT","R2_total_ACT","Active R2 complexes","Active_R2_tot",ngml,visibility,ligandmap);

%% R1, R2 active (surface) 

activeFig5 = VizActiveV121V165P1P2(Routs,"R1_surf_ACT","R1_surf_ACT","R1_surf_ACT","R1_surf_ACT","Active R1 complexes","Active_R1_surf",ngml,visibility,ligandmap);
activeFig6 = VizActiveV121V165(Routs,"R2_surf_ACT","R2_surf_ACT","Active R2 complexes","Active_R2_surf",ngml,visibility,ligandmap);

%% R1, R2 active (internal) 

activeFig7 = VizActiveV121V165P1P2(Routs,"R1_int_ACT","R1_int_ACT","R1_int_ACT","R1_int_ACT","Active R1 complexes","Active_R1_int",ngml,visibility,ligandmap);
activeFig8 = VizActiveV121V165(Routs,"R2_int_ACT","R2_int_ACT","Active R2 complexes","Active_R2_int",ngml,visibility,ligandmap);

%% Ligand Binding on the surface

BndLigFig = VizActiveV121V165P1P2(Routs,"BndV121_surf","BndV165_surf","BndPGF1_surf","BndPGF2_surf","Bound Ligand","BndLigandSurf",ngml,visibility,ligandmap);

%% Total amounts of receptor

TotalRecFig1 = VizActiveV121V165P1P2(Routs,"R1_total","R1_total","R1_total","R1_total","Receptor level","RecLevelTot_R1",ngml,visibility,ligandmap);
TotalRecFig2 = VizActiveV121V165P1P2(Routs,"R1_surf","R1_surf","R1_surf","R1_surf","Receptor level","RecLevelSurf_R1",ngml,visibility,ligandmap);
TotalRecFig3 = VizActiveV121V165P1P2(Routs,"R1_int","R1_int","R1_int","R1_int","Receptor level","RecLevelInt_R1",ngml,visibility,ligandmap);

TotalRecFig4 = VizActiveV121V165P1P2(Routs,"R2_total","R2_total","R2_total","R2_total","Receptor level","RecLevelTot_R2",ngml,visibility,ligandmap);
TotalRecFig5 = VizActiveV121V165P1P2(Routs,"R2_surf","R2_surf","R2_surf","R2_surf","Receptor level","RecLevelSurf_R2",ngml,visibility,ligandmap);
TotalRecFig6 = VizActiveV121V165P1P2(Routs,"R2_int","R2_int","R2_int","R2_int","Receptor level","RecLevelInt_R2",ngml,visibility,ligandmap);

TotalRecFig7 = VizActiveV121V165P1P2(Routs,"N1_total","N1_total","N1_total","N1_total","Receptor level","RecLevelTot_N1",ngml,visibility,ligandmap);
TotalRecFig8 = VizActiveV121V165P1P2(Routs,"N1_surf","N1_surf","N1_surf","N1_surf","Receptor level","RecLevelSurf_N1",ngml,visibility,ligandmap);
TotalRecFig9 = VizActiveV121V165P1P2(Routs,"N1_int","N1_int","N1_int","N1_int","Receptor level","RecLevelInt_N1",ngml,visibility,ligandmap);


%% Fluxes

FluxFig1 = VizFluxes(R1rates,"R1","Flux (#/sec)",ngml,visibility);
FluxFig2 = VizFluxes(R2rates,"R2","Flux (#/sec)",ngml,visibility);
FluxFig3 = VizFluxes(N1rates,"N1","Flux (#/sec)",ngml,visibility);
FluxFig4 = VizFluxes(R1ratesnorm,"R1norm","Flux (relative to no ligand)",ngml,visibility);
FluxFig5 = VizFluxes(R2ratesnorm,"R2norm","Flux (relative to no ligand)",ngml,visibility);
FluxFig6 = VizFluxes(N1ratesnorm,"N1norm","Flux (relative to no ligand)",ngml,visibility);


%% Net Fluxes

NetFluxFig1 = VizNetFluxes(R1rates,"R1",ngml,visibility);
NetFluxFig2 = VizNetFluxes(R2rates,"R2",ngml,visibility);
NetFluxFig3 = VizNetFluxes(N1rates,"N1",ngml,visibility);


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
