close all;
clear all;

model = "All_Lig_Traff_model_Bionetgen_20231016"; %"TraffPhosph"
parameters = base_parameters(model);

% HIGHER VEGFR2 INTERNALIZATION
parameters.kVegfR2Rab5a = parameters.kR2Rab5a * 3;
parameters.kVegfR2N1Rab5a = parameters.kVegfR2Rab5a;
parameters0 = parameters;

% ngml = [0 0.01 0.025 0.05 0.1 0.25 0.5 1 2.5 5 10 25 50 100 150];
ngml = [0 logspace(-2,3,32 )];
minngml = min(ngml(2:end));
maxngml = max(ngml(2:end));
timepts = 0:60:(240*60); %[0 15 30 60 60 240];

ReFitProductionRates = 0 ; % 1 for yes, 0 for no

visibility = 'off';
longfigs = 0; 

%% Surface area of cellular compartments
SA_surface = 1000; % (um2/cell)
SA_rab4    = 950; % (um2/cell)
SA_rab11   = 325; % (um2/cell)

vol_surface = 1e7; % (fL/cell)
vol_rab4    = 11.25; % (fL/cell)
vol_rab11   = 3.75; % (fL/cell)

%% DEFINE LIGAND DOSE & SCATCHARD SCENARIOS TO RUN

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

% runCase = "onlyR1";
% runCase = "onlyR2";
% runCase = "onlyN1";
% runCase = "R1&R2";
% runCase = "R1&N1";
% runCase = "N1&R2";
% runCase = "All3R";

% DimType = "121"; % assuming 1:1 binding of ligands to receptors
% DimType = "LID"; % assuming no pre-dimerization of receptors
% DimType = "DPD"; % assuming only pre-dimerization of receptors and no ligand-induced dimerization
                   % Note: DPD mode is not quite right (because it still permits ligands binding to R monomers without productive further binding)
% DimType = "All"; % all dimerization modes

ScatchardMode = "on"; % when "on", assumes no production, no internalization, only binding (4dC)

scatchCases = ["onlyR1","onlyR2","onlyN1","R1&R2","R1&N1","N1&R2","All3R"];
% scatchCases = ["onlyR1"];
dimCases = ["121","LID","All"];

for j = 1:length(scatchCases)
    runCase = scatchCases(j);

    for k = 1:length(dimCases)
        DimType = dimCases(k);
        
        parameters = parameters0;
        R1surf_target = 1800;
        R2surf_target = 4900; 
        N1surf_target = 68000;

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

%% ACTUALLY RUN THE CASES

disp('running ligand cases')

offset = ((j-1)*length(dimCases)+(k-1))*4*length(ngml);
for i = 1:length(v121a)
    Routs(offset+i) = RunLigandCase(model,speciesInit,parameters,v121a(i),v165a(i),plgf1(i),plgf2(i),timepts);
    disp(strcat('Ligand cases =',num2str((offset+i)*100/(length(scatchCases)*length(dimCases)*length(v121a))),'% complete'))
    [R1rates(offset+i),R2rates(offset+i),N1rates(offset+i)] = CalcRatesL(Routs(offset+i),parameters); % perturbs?
end

    end
end
%% end of calculation loop

outnames = fieldnames(Routs);
ratenames = fieldnames(R1rates);

mapr = [linspace(0.33,1,101)  linspace(.99,0,100)];
mapg = [linspace(0,1,101)   linspace(.99,.33,100)];
mapb = [linspace(0.33,1,101) linspace(0.99,0,100)];
map =[mapr' mapg' mapb'];

% blue,red,green maps for R1 lines
map1 = [linspace(.75,0,length(ngml))];
map2 = [linspace(1,1,length(ngml))];
mapblue =[map1' map1' map2'];
mapred =[map2' map1' map1'];
mapgreen =[map1' map2' map1'];
mapall = [mapblue;mapred;mapgreen];

%% OUTPUT TIMECOURSES FOR THE METRICS

if longfigs == 1
    
disp('creating figures')

% ADD ADDITIONAL FIGURES HERE AS NEEDED

end

% Scatchard plots - VEGFR1

for j = 1:length(scatchCases)  % scatchCases = ["onlyR1","onlyR2","onlyN1"];

a1 = []; a2 = []; a3 = []; 
b1 = []; b2 = []; b3 = []; 
c1 = []; c2 = []; c3 = []; 
d1 = []; d2 = []; d3 = [];
e1 = []; e2 = []; e3 = []; 
f1 = []; f2 = []; f3 = []; 
g1 = []; g2 = []; g3 = []; 
h1 = []; h2 = []; h3 = [];

for i = 1:length(ngml)
    k = 1; % dimCases = ["121","LID","All"];
    offset = ((j-1)*length(dimCases)+(k-1))*4*length(ngml);
	a1 = [a1 Routs(offset+i).BndV121_surf(end)];
    b1 = [b1 (Routs(offset+i).BndV121_surf(end))./(Routs(offset+i).V121_surf(end))];
	c1 = [c1 Routs(offset+i+length(ngml)).BndV165_surf(end)];
    d1 = [d1 (Routs(offset+i+length(ngml)).BndV165_surf(end))./(Routs(offset+i+length(ngml)).V165_surf(end))];
	e1 = [e1 Routs(offset+i+2*length(ngml)).BndPGF1_surf(end)];
    f1 = [f1 (Routs(offset+i+2*length(ngml)).BndPGF1_surf(end))./(Routs(offset+i+2*length(ngml)).PGF1_surf(end))];
	g1 = [g1 Routs(offset+i+3*length(ngml)).BndPGF2_surf(end)];
    h1 = [h1 (Routs(offset+i+3*length(ngml)).BndPGF2_surf(end))./(Routs(offset+i+3*length(ngml)).PGF2_surf(end))];
    k = 2; % dimCases = ["121","LID","All"];
    offset = ((j-1)*length(dimCases)+(k-1))*4*length(ngml);
	a2 = [a2 Routs(offset+i).BndV121_surf(end)];
    b2 = [b2 (Routs(offset+i).BndV121_surf(end))./(Routs(offset+i).V121_surf(end))];
	c2 = [c2 Routs(offset+i+length(ngml)).BndV165_surf(end)];
    d2 = [d2 (Routs(offset+i+length(ngml)).BndV165_surf(end))./(Routs(offset+i+length(ngml)).V165_surf(end))];
	e2 = [e2 Routs(offset+i+2*length(ngml)).BndPGF1_surf(end)];
    f2 = [f2 (Routs(offset+i+2*length(ngml)).BndPGF1_surf(end))./(Routs(offset+i+2*length(ngml)).PGF1_surf(end))];
	g2 = [g2 Routs(offset+i+3*length(ngml)).BndPGF2_surf(end)];
    h2 = [h2 (Routs(offset+i+3*length(ngml)).BndPGF2_surf(end))./(Routs(offset+i+3*length(ngml)).PGF2_surf(end))];
    k = 3; % dimCases = ["121","LID","All"];
    offset = ((j-1)*length(dimCases)+(k-1))*4*length(ngml);
	a3 = [a3 Routs(offset+i).BndV121_surf(end)];
    b3 = [b3 (Routs(offset+i).BndV121_surf(end))./(Routs(offset+i).V121_surf(end))];
	c3 = [c3 Routs(offset+i+length(ngml)).BndV165_surf(end)];
    d3 = [d3 (Routs(offset+i+length(ngml)).BndV165_surf(end))./(Routs(offset+i+length(ngml)).V165_surf(end))];
	e3 = [e3 Routs(offset+i+2*length(ngml)).BndPGF1_surf(end)];
    f3 = [f3 (Routs(offset+i+2*length(ngml)).BndPGF1_surf(end))./(Routs(offset+i+2*length(ngml)).PGF1_surf(end))];
	g3 = [g3 Routs(offset+i+3*length(ngml)).BndPGF2_surf(end)];
    h3 = [h3 (Routs(offset+i+3*length(ngml)).BndPGF2_surf(end))./(Routs(offset+i+3*length(ngml)).PGF2_surf(end))];
end

figname = strcat("ScatchardFig",num2str(j));
eval(strcat(figname," = figure('visible',visibility);"));
ax1 = subplot(1,4,1);
plot(a1,b1,'LineWidth',2); hold on;
plot(a2,b2,'LineWidth',2); 
plot(a3,b3,'LineWidth',2); 
legend(dimCases);
title(ax1, 'VEGF121');
xlabel(ax1,'Bound');
ylabel(ax1,'Bound/Free');
a = [a1' b1' a2' b2' a3' b3'];
csvname = strcat("outputDataDr0/Scatchard_",scatchCases(j),"_a.csv");
csvwrite(csvname,a);
  ax2 = subplot(1,4,2);
  plot(c1,d1,'LineWidth',2); hold on;
  plot(c2,d2,'LineWidth',2); 
  plot(c3,d3,'LineWidth',2); 
  title(ax2, 'VEGF165');
  xlabel(ax2,'Bound');
  ylabel(ax2,'Bound/Free');
  b = [c1' d1' c2' d2' c3' d3'];
  csvname = strcat("outputDataDr0/Scatchard_",scatchCases(j),"_b.csv");
  csvwrite(csvname,b);
    ax3 = subplot(1,4,3);
    plot(e1,f1,'LineWidth',2); hold on;
    plot(e2,f2,'LineWidth',2); 
    plot(e3,f3,'LineWidth',2); 
    title(ax3, 'PlGF1');
    xlabel(ax3,'Bound');
    ylabel(ax3,'Bound/Free');
    c = [e1' f1' e2' f2' e3' f3'];
    csvname = strcat("outputDataDr0/Scatchard_",scatchCases(j),"_c.csv");
    csvwrite(csvname,c);
      ax4 = subplot(1,4,4);
      plot(g1,h1,'LineWidth',2); hold on;
      plot(g2,h2,'LineWidth',2); 
      plot(g3,h3,'LineWidth',2); 
      title(ax4, 'PlGF2');
      xlabel(ax4,'Bound');
      ylabel(ax4,'Bound/Free');
      d = [g1' h1' g2' h2' g3' h3'];
      csvname = strcat("outputDataDr0/Scatchard_",scatchCases(j),"_d.csv");
      csvwrite(csvname,d);
set(eval(figname), 'Position',[0 0 1000 200])
set(eval(figname), 'PaperUnits', 'inches');
set(eval(figname), 'PaperSize', [4 6]);
exportgraphics(eval(figname),strcat("outputFigsDr0/Scatchard_",scatchCases(j),".png"),'Resolution',300);

end




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
