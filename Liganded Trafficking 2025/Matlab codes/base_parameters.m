function parameters = base_parameters(model)

switch model
    case "All_Lig_Traff_model_Bionetgen_20231016"

% set up the parameter structure (using a consistent order)
a = fileread("Paramlist.txt");
pn = strsplit(a);
for i=1:length(pn)
    parameters.(pn{i}) = 0;
end


%% INITIAL CONDITIONS

% Initial Ligand Levels (set to zero for receptor-only pre-simulation)
parameters.VEGF165a_0 = 0; %6843*1000;  % = 50ng/mL  
parameters.PLGF1_0 = 0;    %10140*1000; % = 50ng/mL 
parameters.PLGF2_0 = 0;    %8702*1000;  % = 50ng/mL 
parameters.VEGF121a_0 = 0; %10750*1000;	% = 50ng/mL

% Initial Receptor Levels
% parameter values from Imoukhuede & Popel. Experimental Cell Research (2011) https://doi.org/10.1016/j.yexcr.2010.12.014
% All in molecules per cell
parameters.VEGFR1_0 = 18000; % total R1 calculated from 1800 rec/cell on surface/total fraction of 0.1
parameters.VEGFR2_0 =  9607; % total R2 calculated from 4900 rec/cell on surface/total fraction of 0.51
parameters.NRP1_0   = 91891; % total R1 calculated from 68000 rec/cell on surface/total fraction of 0.74

parameters.M_0 = 0; % set to zero for these simulations;

%% 1. RECEPTOR-RECEPTOR COUPLING/UNCOUPLING PARAMETERS (kon = (molec/cell)-1s-1; koff ~ s-1)

%% Surface area of cellular compartments
SA_surface = 1000; % (um2/cell)
SA_rab4    = 950; % (um2/cell)
SA_rab11   = 325; % (um2/cell)

%% Receptor-Receptor association/dissociation parameters 
% Units: kon = (receptors/cell)-1.s-1; koff = s-1
kR1R1on = 8e-4; % units = (receptors/um2)-1.s-1
parameters.kR1R1on_surf  = kR1R1on/SA_surface;
parameters.kR1R1on_rab4  = kR1R1on/SA_rab4;
parameters.kR1R1on_rab11 = kR1R1on/SA_rab11;
parameters.kR1R1off = 0.01;
kR2R2on = 2e-3; % units = (receptors/um2)-1.s-1
parameters.kR2R2on_surf  = kR2R2on/SA_surface;
parameters.kR2R2on_rab4  = kR2R2on/SA_rab4;
parameters.kR2R2on_rab11 = kR2R2on/SA_rab11;
parameters.kR2R2off = 0.01;

% parameters.kR1R1on_surf = 8e-07;
% parameters.kR1R1on_rab4 = 8.42e-07;
% parameters.kR1R1on_rab11 = 2.46e-06;
% parameters.kR1R1off = 0.01;
% parameters.kR2R2on_surf = 2e-06;
% parameters.kR2R2on_rab4 = 2.11e-06;
% parameters.kR2R2on_rab11 = 6.15e-06;
% parameters.kR2R2off = 0.01;

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

%% 2. LIGAND-RECEPTOR BINDING/UNBINDING PARAMETERS (kon = (molec/um2)-1s-1; koff ~ s-1)

% rate constants for ligands binding to receptors
% parameter values based on estimates from Clegg & Mac Gabhann (2007) PLOS.Comp.Biol https://doi.org/10.1371/journal.pcbi.1005445

% LIGANDS BINDING TO RECEPTOR MONOMERS or DIMERS

%% Volume of cellular compartments
vol_surface = 1e7; % (fL/cell)
vol_rab4    = 11.25; % (fL/cell)
vol_rab11   = 3.75; % (fL/cell)

parameters.kvegfA_R1on_surf = 1.245e-09;
parameters.kvegfA_R1on_rab4 = parameters.kvegfA_R1on_surf * vol_surface/vol_rab4;
parameters.kvegfA_R1on_rab11 = parameters.kvegfA_R1on_surf * vol_surface/vol_rab11;
parameters.kvegfA_R1off = 0.0224;

% parameters.kvegfA_R1on_surf = 1.245e-09;
% parameters.kvegfA_R1on_rab4 = 0.001107;
% parameters.kvegfA_R1on_rab11 = 0.003321;
% parameters.kvegfA_R1off = 0.0224;

parameters.kN1R1_vegfAon_surf = parameters.kvegfA_R1on_surf;
parameters.kN1R1_vegfAon_rab4 = parameters.kvegfA_R1on_rab4;
parameters.kN1R1_vegfAon_rab11 = parameters.kvegfA_R1on_rab11;
parameters.kN1R1_vegfAoff = parameters.kvegfA_R1off;

parameters.kplgf_R1on_surf = 6.227e-11;
parameters.kplgf_R1on_rab4 = parameters.kplgf_R1on_surf * vol_surface/vol_rab4;
parameters.kplgf_R1on_rab11 = parameters.kplgf_R1on_surf * vol_surface/vol_rab11;
parameters.kplgf_R1off = 0.0132;

% parameters.kplgf_R1on_surf = 6.227e-11;
% parameters.kplgf_R1on_rab4 = 5.535e-05;
% parameters.kplgf_R1on_rab11 = 0.0001661;
% parameters.kplgf_R1off = 0.0132;

parameters.kN1R1_plgfon_surf = parameters.kplgf_R1on_surf;
parameters.kN1R1_plgfon_rab4 = parameters.kplgf_R1on_rab4;
parameters.kN1R1_plgfon_rab11 = parameters.kplgf_R1on_rab11;
parameters.kN1R1_plgfoff = parameters.kplgf_R1off;

parameters.kvegfA_R2on_surf = 4.151e-10;
parameters.kvegfA_R2on_rab4 = parameters.kvegfA_R2on_surf * vol_surface/vol_rab4;
parameters.kvegfA_R2on_rab11 = parameters.kvegfA_R2on_surf * vol_surface/vol_rab11;
parameters.kvegfA_R2off = 0.0224;

parameters.kvegfA_N1on_surf = 2.076e-11;
parameters.kvegfA_N1on_rab4 = parameters.kvegfA_N1on_surf * vol_surface/vol_rab4;
parameters.kvegfA_N1on_rab11 = parameters.kvegfA_N1on_surf * vol_surface/vol_rab11;
parameters.kvegfA_N1off = 0.0173;

parameters.kplgf_N1on_surf = 4.151e-13;
parameters.kplgf_N1on_rab4 = parameters.kplgf_N1on_surf * vol_surface/vol_rab4;
parameters.kplgf_N1on_rab11 = parameters.kplgf_N1on_surf * vol_surface/vol_rab11;
parameters.kplgf_N1off = 0.0224;

% parameters.kvegfA_R2on_surf = 4.151e-10;
% parameters.kvegfA_R2on_rab4 = 0.000369;
% parameters.kvegfA_R2on_rab11 = 0.001107;
% parameters.kvegfA_R2off = 0.0224;
% parameters.kvegfA_N1on_surf = 2.076e-11;
% parameters.kvegfA_N1on_rab4 = 1.845e-05;
% parameters.kvegfA_N1on_rab11 = 5.535e-05;
% parameters.kvegfA_N1off = 0.0173;
% parameters.kplgf_N1on_surf = 4.151e-13;
% parameters.kplgf_N1on_rab4 = 3.69e-07;
% parameters.kplgf_N1on_rab11 = 1.107e-06;
% parameters.kplgf_N1off = 0.0224;

parameters.kR1vegfA_R1on_surf = 5.307e-04 ;
parameters.kR1vegfA_R1on_rab4 = parameters.kR1vegfA_R1on_surf * SA_surface/SA_rab4;
parameters.kR1vegfA_R1on_rab11 = parameters.kR1vegfA_R1on_surf * SA_surface/SA_rab11;
parameters.kR1plgf_R1on_surf = 5.409e-04 ;
parameters.kR1plgf_R1on_rab4 = parameters.kR1plgf_R1on_surf * SA_surface/SA_rab4;
parameters.kR1plgf_R1on_rab11 = parameters.kR1plgf_R1on_surf * SA_surface/SA_rab11;
parameters.kR2vegfA_R2on_surf = 1.911e-04 ;
parameters.kR2vegfA_R2on_rab4 = parameters.kR2vegfA_R2on_surf * SA_surface/SA_rab4;
parameters.kR2vegfA_R2on_rab11 = parameters.kR2vegfA_R2on_surf * SA_surface/SA_rab11;

% parameters.kR1vegfA_R1on_surf = 5.307e-04;
% parameters.kR1vegfA_R1on_rab4 = 5.586e-04;
% parameters.kR1vegfA_R1on_rab11 = 0.001633;
% parameters.kR1plgf_R1on_surf = 5.409e-04;
% parameters.kR1plgf_R1on_rab4 = 5.693e-04;
% parameters.kR1plgf_R1on_rab11 = 0.001664;
% parameters.kR2vegfA_R2on_surf = 1.911e-04;
% parameters.kR2vegfA_R2on_rab4 = 2.011e-04;
% parameters.kR2vegfA_R2on_rab11 = 5.879e-04;

parameters.kN1vegfA_R2on_surf = 2.5e-6 ; %/10; % parameters.kR2vegfA_R2on_surf;
parameters.kN1vegfA_R2on_rab4 = parameters.kN1vegfA_R2on_surf * SA_surface/SA_rab4; %parameters.kR2vegfA_R2on_rab4;
parameters.kN1vegfA_R2on_rab11 = parameters.kN1vegfA_R2on_surf * SA_surface/SA_rab4; %parameters.kR2vegfA_R2on_rab11;
parameters.kN1vegfA_R2off = parameters.kvegfA_R2off;

parameters.kR2vegfA_N1on_surf = 2.5e-6 ; %/10; % 1.420e-05 (parameters.kN1vegfA_N1on_surf);
parameters.kR2vegfA_N1on_rab4 = parameters.kR2vegfA_N1on_surf * SA_surface/SA_rab4;
parameters.kR2vegfA_N1on_rab11 = parameters.kR2vegfA_N1on_surf * SA_surface/SA_rab11;
parameters.kR2vegfA_N1off = parameters.kvegfA_N1off;
% parameters.kR2vegfA_N1on_surf = 1.420e-05;
% parameters.kR2vegfA_N1on_rab4 = 1.494e-05;
% parameters.kR2vegfA_N1on_rab11 = 4.368e-05;
% parameters.kR2vegfA_N1off = 0.0173;
parameters.kN1vegfA_N1on_surf = 1.420e-05; %
parameters.kN1vegfA_N1on_rab4 = parameters.kN1vegfA_N1on_surf * SA_surface/SA_rab4;
parameters.kN1vegfA_N1on_rab11 = parameters.kN1vegfA_N1on_surf * SA_surface/SA_rab11;
parameters.kN1vegfA_N1off = parameters.kvegfA_N1off;

parameters.kN1plgf2_N1on_surf = 1.405e-05;
parameters.kN1plgf2_N1on_rab4 = parameters.kN1plgf2_N1on_surf * SA_surface/SA_rab4;
parameters.kN1plgf2_N1on_rab11 = parameters.kN1plgf2_N1on_surf * SA_surface/SA_rab11;
parameters.kN1plgf2_N1off = parameters.kplgf_N1off;
% parameters.kN1plgf2_N1on_surf = 1.405e-05;
% parameters.kN1plgf2_N1on_rab4 = 1.479e-05;
% parameters.kN1plgf2_N1on_rab11 = 4.323e-05;
% parameters.kN1plgf2_N1off = 0.0173;

parameters.kdeltaR1R1_V = 0.00144;
parameters.kdeltaR1R1_P = 0.00144;
parameters.kdeltaR2R2 = 0.01;
parameters.kdeltaVEGFA_R = 0.95528;
parameters.kdeltaPLGF_R = 0.97354;

parameters.kvegfAMon = 2.658e-08;
parameters.kvegfAMoff = 0.01;
parameters.kplgfMon = 3.65e-08;
parameters.kplgfMoff = 0.001;
parameters.kMplgf_R1on = 3.987e-07;
parameters.kMplgf_R1off = 0.00035;
parameters.kMvegfA_R1on = 4.98e-06;
parameters.kMvegfA_R1off = 0.001;
parameters.kMvegfA_R2on = 1.66e-06;
parameters.kMvegfA_R2off = 0.001;
parameters.kR1vegfA_Mon = 2.658e-08;
parameters.kR1vegfA_Moff = 0.01;
parameters.kR1plgf_Mon = 3.65e-08;
parameters.kR1plgf_Moff = 0.001;
parameters.kR2vegfA_Mon = 2.658e-08;
parameters.kR2vegfA_Moff = 0.01;

%% 4. TRAFFICKING PARAMETERS 
% (all first order processes, all s-1)

%% INTERNALIZATION
 
parameters.kR1Rab5a = 1.347E-2;
parameters.kR2Rab5a = 2.304E-4;
parameters.kN1Rab5a = 2.685E-4;
 
parameters.kR1N1Rab5a = parameters.kR1Rab5a;
parameters.kR2N1Rab5a = parameters.kR2Rab5a;
 
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

%% RECYCLING via Rab4
 
parameters.kR1Rab4a = 5.367E-4;
parameters.kR2Rab4a = 1.195E-6;
parameters.kN1Rab4a = 2.127E-2;
 
parameters.kR1N1Rab4a = parameters.kR1Rab4a;
parameters.kR2N1Rab4a = parameters.kR2Rab4a;
 
parameters.kVegfR1Rab4a = parameters.kR1Rab4a;
parameters.kVegfR1N1Rab4a = parameters.kVegfR1Rab4a;
parameters.kPlgfR1Rab4a = parameters.kR1Rab4a;
parameters.kPlgfR1N1Rab4a = parameters.kPlgfR1Rab4a;
parameters.kVegfR2Rab4a = parameters.kR2Rab4a;
parameters.kVegfR2N1Rab4a = parameters.kVegfR2Rab4a;
 
parameters.kVegfN1Rab4a = parameters.kN1Rab4a;
parameters.kPlgfN1Rab4a = parameters.kN1Rab4a;
 
parameters.kVegfRab4a = 0;
parameters.kPlgfRab4a = 0;
%% TRANSFER Rab4=>Rab11
 
parameters.kR1Rab4at11a = 5.909E-4;
parameters.kR2Rab4at11a = 1.494E-6;
parameters.kN1Rab4at11a = 7.025E-2;
 
parameters.kR1N1Rab4at11a = parameters.kR1Rab4at11a;
parameters.kR2N1Rab4at11a = parameters.kR2Rab4at11a;
 
parameters.kVegfR1Rab4at11a = parameters.kR1Rab4at11a;
parameters.kVegfR1N1Rab4at11a = parameters.kVegfR1Rab4at11a;
parameters.kPlgfR1Rab4at11a = parameters.kR1Rab4at11a;
parameters.kPlgfR1N1Rab4at11a = parameters.kPlgfR1Rab4at11a;
parameters.kVegfR2Rab4at11a = parameters.kR2Rab4at11a;
parameters.kVegfR2N1Rab4at11a = parameters.kVegfR2Rab4at11a;
 
parameters.kVegfN1Rab4at11a = parameters.kN1Rab4at11a;
parameters.kPlgfN1Rab4at11a = parameters.kN1Rab4at11a;
 
parameters.kVegfRab4at11a = 0;
parameters.kPlgfRab4at11a = 0;
 
%% RECYCLING via Rab11
 
parameters.kR1Rab11a = 9.999E-2;
parameters.kR2Rab11a = 8.933E-2;
parameters.kN1Rab11a = 7.856E-4;
 
parameters.kR1N1Rab11a = parameters.kR1Rab11a;
parameters.kR2N1Rab11a = parameters.kR2Rab11a;

parameters.kVegfR1Rab11a = parameters.kR1Rab11a;
parameters.kVegfR1N1Rab11a = parameters.kVegfR1Rab11a;
parameters.kPlgfR1Rab11a = parameters.kR1Rab11a;
parameters.kPlgfR1N1Rab11a = parameters.kPlgfR1Rab11a;
parameters.kVegfR2Rab11a = parameters.kR2Rab11a;
parameters.kVegfR2N1Rab11a = parameters.kVegfR2Rab11a;
 
parameters.kVegfN1Rab11a = parameters.kN1Rab11a;
parameters.kPlgfN1Rab11a = parameters.kN1Rab11a;
 
parameters.kVegfRab11a = 0;
parameters.kPlgfRab11a = 0;

%% DEGRADATION
 
parameters.kR1Rab4at7a = 2.293E-4;
parameters.kR2Rab4at7a = 2.337E-4;
parameters.kN1Rab4at7a = 1.228E-6;
 
parameters.kR1N1Rab4at7a = parameters.kR1Rab4at7a;
parameters.kR2N1Rab4at7a = parameters.kR2Rab4at7a;

parameters.kVegfR1Rab4at7a = parameters.kR1Rab4at7a;
parameters.kVegfR1N1Rab4at7a = parameters.kVegfR1Rab4at7a;
parameters.kPlgfR1Rab4at7a = parameters.kR1Rab4at7a;
parameters.kPlgfR1N1Rab4at7a = parameters.kPlgfR1Rab4at7a;
parameters.kVegfR2Rab4at7a = parameters.kR2Rab4at7a;
parameters.kVegfR2N1Rab4at7a = parameters.kVegfR2Rab4at7a;
 
parameters.kVegfN1Rab4at7a = parameters.kN1Rab4at7a;
parameters.kPlgfN1Rab4at7a = parameters.kN1Rab4at7a;
 
%% LIGAND DEGRADATION
 
parameters.kVegfRab4at7a = 0.012;
parameters.kPlgfRab4at7a = 0.012;
 
 
%% RECEPTOR PRODUCTION
parameters.kR1prod = 4.101;
parameters.kR2prod = 1.114;
parameters.kN1prod = 0.459;

    

end
end
