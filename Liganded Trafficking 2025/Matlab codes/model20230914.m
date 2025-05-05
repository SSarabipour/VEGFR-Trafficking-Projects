function [R1_total, R2_total, N1_total, R1_surf, R2_surf, N1_surf, R1_int, R2_int, N1_int, R1_Rab4, R2_Rab4, N1_Rab4, R1_Rab11, R2_Rab11, N1_Rab11, ...
    R1R1_165a_total, R2R2_165a_total, R1R1_165a_surf, R2R2_165a_surf, R1R1_165a_int, R2R2_165a_int, R1R1_plgf1_total, R1R1_plgf1_surf, R1R1_plgf1_int...
    R1R1_165a_total_ACT, R2R2_165a_total_ACT, R1R1_165a_surf_ACT, R2R2_165a_surf_ACT, R1R1_165a_int_ACT, R2R2_165a_int_ACT, R1R1_plgf1_total_ACT, R1R1_plgf1_surf_ACT, R1R1_plgf1_int_ACT, ...
    R1R1_121a_total_ACT, R2R2_121a_total_ACT, R1R1_121a_surf_ACT, R2R2_121a_surf_ACT, R1R1_121a_int_ACT, R2R2_121a_int_ACT, R1R1_plgf2_total_ACT, R1R1_plgf2_surf_ACT, R1R1_plgf2_int_ACT, ...
	V165_surf,V165_int,V165_Rab4,V165_Rab11,V121_surf,V121_int,V121_Rab4,V121_Rab11, ...
    PGF1_surf,PGF1_int,PGF1_Rab4,PGF1_Rab11,PGF2_surf,PGF2_int,PGF2_Rab4,PGF2_Rab11, ligand_surf, ligand_int,...
    R1_empty_surf,R2_empty_surf,N1_empty_surf,R1_empty_Rab4,R2_empty_Rab4,N1_empty_Rab4,R1_empty_Rab11,R2_empty_Rab11,N1_empty_Rab11] = model20230914(runModel,inits,parameters,vgf165,vgf121,pgf1,pgf2,timepts) 

% This code runs one simulation of the model
% OUTPUTS: R1_total, R2_total, N1_total = whole-cell total receptor levels (i.e. ligand and unliganded)
%         R1_surf, R2_surf, N1_surf = surface total receptor levels (i.e. ligand and unliganded)
%         R1_surf, R2_surf, N1_surf = surface total receptor levels (i.e. ligand and unliganded)
%         R1R1_165a_surf,  R2R2_165a_surf = surface total active receptors
%         R1R1_165a_total, R2R2_165a_total = whole cell total active receptors
% INPUTS: vgf, pgf = vegf and plgf concentrations for these ligand
%         inits = initial conditions
%         parameters = parameter values
%         timepts = which timepoint results to return to the main code

%% INITIAL CONCENTRATIONS OF LIGANDS

inits(1) = inits(1) + vgf165;
inits(2) = inits(2) + vgf121;
inits(3) = inits(3) + pgf1;
inits(4) = inits(4) + pgf2;

% inits(1) = inits(1) + vgf;
% inits(4) = inits(4) + vgf;
% initial4(2) = initial4(2) + parameters.mVEGF165b_0; %monomeric VEGF-A165b
% initial4(3) = initial4(3) + parameters.dVEGF165b_0;  %dimeric VEGF-A165b
% inits(2) = inits(2) + pgf;
% inits(3) = inits(3) + pgf;
%inits(3) = inits(3) + pgf2;
%inits(4) = inits(4) + vgf121a;
% initial4(6) = initial4(6) + parameters.PLGF2_0;

%   parameters.VEGF165a_0 = vgf; %50ng/mL % (molecules/um^2)
%   parameters.PLGF1_0 = pgf;

%% RUN MODEL
TimeLen = 3600*4;
TimeStep = 6;
Timepoints  = 0:TimeStep:TimeLen; 
% [timepoints, species_out, observables_out] = Ligand_Traff_model_14Sep2023(timeLen, timestp_sim, parameters, inits);
[~, timepoints, species_out, observables_out] = eval(append(runModel,"(Timepoints', inits, parameters, 1)"));
ConcOut=CalcOutputsL(observables_out);

%% OUTPUTS
for i = 1:length(timepts)
    j=timepts(i)*10+1; % every six seconds, so number of minutes times ten (plus one for time zero)
    R1_total(i) = ConcOut.R1_total(j) ;
    R2_total(i) = ConcOut.R2_total(j) ;
    N1_total(i) = ConcOut.N1_total(j) ;
    R1_surf(i)  = ConcOut.R1_surf(j) ;
    R2_surf(i)  = ConcOut.R2_surf(j) ;
    N1_surf(i)  = ConcOut.N1_surf(j) ;
    R1_int(i)  = ConcOut.R1_int(j) ;
    R2_int(i)  = ConcOut.R2_int(j) ;
    N1_int(i)  = ConcOut.N1_int(j) ;
    R1_Rab4(i)  = ConcOut.R1_Rab4(j) ;
    R2_Rab4(i)  = ConcOut.R2_Rab4(j) ;
    N1_Rab4(i)  = ConcOut.N1_Rab4(j) ;
    R1_Rab11(i)  = ConcOut.R1_Rab11(j) ;
    R2_Rab11(i)  = ConcOut.R2_Rab11(j) ;
    N1_Rab11(i)  = ConcOut.N1_Rab11(j) ;
    
    R1R1_165a_total(i) = ConcOut.R1R1_165a_total(j) ;
    R2R2_165a_total(i) = ConcOut.R2R2_165a_total(j) ;
    R1R1_165a_surf(i)  = ConcOut.R1R1_165a_surf(j) ;
    R2R2_165a_surf(i)  = ConcOut.R2R2_165a_surf(j) ;
    R1R1_165a_int(i)  = ConcOut.R1R1_165a_int(j) ;
    R2R2_165a_int(i)  = ConcOut.R2R2_165a_int(j) ;
    R1R1_plgf1_total(i) = ConcOut.R1R1_plgf1_total(j) ;
    R1R1_plgf1_surf(i)  = ConcOut.R1R1_plgf1_surf(j) ;
    R1R1_plgf1_int(i)  = ConcOut.R1R1_plgf1_int(j) ;

    R1R1_165a_total_ACT(i) = ConcOut.R1R1_165a_total_ACT(j) ;
    R2R2_165a_total_ACT(i) = ConcOut.R2R2_165a_total_ACT(j) ;
    R1R1_165a_surf_ACT(i)  = ConcOut.R1R1_165a_surf_ACT(j) ;
    R2R2_165a_surf_ACT(i)  = ConcOut.R2R2_165a_surf_ACT(j) ;
    R1R1_165a_int_ACT(i)  = ConcOut.R1R1_165a_int_ACT(j) ;
    R2R2_165a_int_ACT(i)  = ConcOut.R2R2_165a_int_ACT(j) ;
    R1R1_plgf1_total_ACT(i) = ConcOut.R1R1_plgf1_total_ACT(j) ;
    R1R1_plgf1_surf_ACT(i)  = ConcOut.R1R1_plgf1_surf_ACT(j) ;
    R1R1_plgf1_int_ACT(i)  = ConcOut.R1R1_plgf1_int_ACT(j) ;

	R1R1_121a_total_ACT(i) = ConcOut.R1R1_121a_total_ACT(j) ;
    R2R2_121a_total_ACT(i) = ConcOut.R2R2_121a_total_ACT(j) ;
    R1R1_121a_surf_ACT(i)  = ConcOut.R1R1_121a_surf_ACT(j) ;
    R2R2_121a_surf_ACT(i)  = ConcOut.R2R2_121a_surf_ACT(j) ;
    R1R1_121a_int_ACT(i)  = ConcOut.R1R1_121a_int_ACT(j) ;
    R2R2_121a_int_ACT(i)  = ConcOut.R2R2_121a_int_ACT(j) ;
    R1R1_plgf2_total_ACT(i) = ConcOut.R1R1_plgf2_total_ACT(j) ;
    R1R1_plgf2_surf_ACT(i)  = ConcOut.R1R1_plgf2_surf_ACT(j) ;
    R1R1_plgf2_int_ACT(i)  = ConcOut.R1R1_plgf2_int_ACT(j) ;

    numpercell_to_pM_surf  = (1.0e4)/(6.022)/1.0e7; % => (#/cell) * (1e12 pmol/mol) * (1e15 fL/L) /(Na #/mol) / (fL/cell) 
    numpercell_to_pM_int   = (1.0e4)/(6.022)/15; 
    numpercell_to_pM_rab4  = (1.0e4)/(6.022)/11.25; 
    numpercell_to_pM_rab11 = (1.0e4)/(6.022)/3.75; 

    V165_surf(i) = ConcOut.V165_surf(j) * numpercell_to_pM_surf;
    V165_int(i)  = ConcOut.V165_int(j) * numpercell_to_pM_int;
    V165_Rab4(i) = ConcOut.V165_Rab4(j) * numpercell_to_pM_rab4;
    V165_Rab11(i)= ConcOut.V165_Rab11(j) * numpercell_to_pM_rab11;
    V121_surf(i) = ConcOut.V121_surf(j) * numpercell_to_pM_surf;
    V121_int(i)  = ConcOut.V121_int(j) * numpercell_to_pM_int;
    V121_Rab4(i) = ConcOut.V121_Rab4(j) * numpercell_to_pM_rab4;
    V121_Rab11(i)= ConcOut.V121_Rab11(j) * numpercell_to_pM_rab11;
    PGF1_surf(i) = ConcOut.PGF1_surf(j) * numpercell_to_pM_surf;
    PGF1_int(i)  = ConcOut.PGF1_int(j) * numpercell_to_pM_int;
    PGF1_Rab4(i) = ConcOut.PGF1_Rab4(j) * numpercell_to_pM_rab4;
    PGF1_Rab11(i)= ConcOut.PGF1_Rab11(j) * numpercell_to_pM_rab11;
    PGF2_surf(i) = ConcOut.PGF2_surf(j) * numpercell_to_pM_surf;
    PGF2_int(i)  = ConcOut.PGF2_int(j) * numpercell_to_pM_int;
    PGF2_Rab4(i) = ConcOut.PGF2_Rab4(j) * numpercell_to_pM_rab4;
    PGF2_Rab11(i)= ConcOut.PGF2_Rab11(j) * numpercell_to_pM_rab11;
    
    ligand_surf(i) = ConcOut.ligand_surf(j) * numpercell_to_pM_surf;
    ligand_int(i) = ConcOut.ligand_int(j) * numpercell_to_pM_int;

    R1_empty_surf(i) = ConcOut.R1_empty_surf(j);
    R2_empty_surf(i) = ConcOut.R2_empty_surf(j);
    N1_empty_surf(i) = ConcOut.N1_empty_surf(j);
    R1_empty_Rab4(i) = ConcOut.R1_empty_Rab4(j);
    R2_empty_Rab4(i) = ConcOut.R2_empty_Rab4(j);
    N1_empty_Rab4(i) = ConcOut.N1_empty_Rab4(j);
    R1_empty_Rab11(i) = ConcOut.R1_empty_Rab11(j);
    R2_empty_Rab11(i) = ConcOut.R2_empty_Rab11(j);
    N1_empty_Rab11(i) = ConcOut.N1_empty_Rab11(j);
end

end