function ConcOut=CalcOutputsL(observables_out);

ConcOut.R1_surf = observables_out(:,1)';
ConcOut.R2_surf = observables_out(:,2)';
ConcOut.N1_surf = observables_out(:,3)';

ConcOut.R1_Rab4 = observables_out(:,4)';
ConcOut.R2_Rab4 = observables_out(:,5)';
ConcOut.N1_Rab4 = observables_out(:,6)';

ConcOut.R1_Rab11= observables_out(:,7)';
ConcOut.R2_Rab11= observables_out(:,8)';
ConcOut.N1_Rab11= observables_out(:,9)';

%ConcOut.R1_Rab7= observables_out(:,)'; % Rab7 here needs to be renamed; should be Golgi
%ConcOut.R2_Rab7= observables_out(:,10)';% Rab7 here needs to be renamed; should be Golgi
%ConcOut.N1_Rab7= observables_out(:,33)';% Rab7 here needs to be renamed; should be Golgi

ConcOut.R1_total = ConcOut.R1_surf+ConcOut.R1_Rab4+ConcOut.R1_Rab11;% Rab7 here needs to be renamed; should be Golgi
ConcOut.R2_total = ConcOut.R2_surf+ConcOut.R2_Rab4+ConcOut.R2_Rab11;% Rab7 here needs to be renamed; should be Golgi
ConcOut.N1_total = ConcOut.N1_surf+ConcOut.N1_Rab4+ConcOut.N1_Rab11;% Rab7 here needs to be renamed; should be Golgi

ConcOut.R1_int   =                 ConcOut.R1_Rab4+ConcOut.R1_Rab11;% Rab7 here needs to be renamed; should be Golgi
ConcOut.R2_int   =                 ConcOut.R2_Rab4+ConcOut.R2_Rab11;% Rab7 here needs to be renamed; should be Golgi
ConcOut.N1_int   =                 ConcOut.N1_Rab4+ConcOut.N1_Rab11;% Rab7 here needs to be renamed; should be Golgi


ConcOut.R2R2_121a_surf = observables_out(:,25)';
ConcOut.R2R2_165a_surf = observables_out(:,26)';

ConcOut.R1R1_121a_surf = observables_out(:,27)';
ConcOut.R1R1_165a_surf = observables_out(:,28)';

ConcOut.R1R1_plgf1_surf = observables_out(:,29)';
ConcOut.R1R1_plgf2_surf = observables_out(:,30)';



ConcOut.R2R2_121a_Rab4 = observables_out(:,37)';
ConcOut.R2R2_165a_Rab4 = observables_out(:,38)';

ConcOut.R1R1_121a_Rab4 = observables_out(:,39)';
ConcOut.R1R1_165a_Rab4 = observables_out(:,40)';

ConcOut.R1R1_plgf1_Rab4 = observables_out(:,41)';
ConcOut.R1R1_plgf2_Rab4 = observables_out(:,42)';



ConcOut.R2R2_121a_Rab11 = observables_out(:,49)';
ConcOut.R2R2_165a_Rab11 = observables_out(:,50)';

ConcOut.R1R1_121a_Rab11 = observables_out(:,51)';
ConcOut.R1R1_165a_Rab11 = observables_out(:,52)';

ConcOut.R1R1_plgf1_Rab11 = observables_out(:,53)';
ConcOut.R1R1_plgf2_Rab11 = observables_out(:,54)';

ConcOut.R2R2_121a_total = ConcOut.R2R2_121a_surf+ConcOut.R2R2_121a_Rab4+ConcOut.R2R2_121a_Rab11;
ConcOut.R2R2_165a_total = ConcOut.R2R2_165a_surf+ConcOut.R2R2_165a_Rab4+ConcOut.R2R2_165a_Rab11;

ConcOut.R1R1_121a_total = ConcOut.R1R1_121a_surf+ConcOut.R1R1_121a_Rab4+ConcOut.R1R1_121a_Rab11;
ConcOut.R1R1_165a_total = ConcOut.R1R1_165a_surf+ConcOut.R1R1_165a_Rab4+ConcOut.R1R1_165a_Rab11;

ConcOut.R1R1_plgf1_total = ConcOut.R1R1_plgf1_surf+ConcOut.R1R1_plgf1_Rab4+ConcOut.R1R1_plgf1_Rab11;
ConcOut.R1R1_plgf2_total = ConcOut.R1R1_plgf2_surf+ConcOut.R1R1_plgf2_Rab4+ConcOut.R1R1_plgf2_Rab11;

ConcOut.R2R2_121a_int = ConcOut.R2R2_121a_Rab4+ConcOut.R2R2_121a_Rab11;
ConcOut.R2R2_165a_int = ConcOut.R2R2_165a_Rab4+ConcOut.R2R2_165a_Rab11;

ConcOut.R1R1_121a_int = ConcOut.R1R1_121a_Rab4+ConcOut.R1R1_121a_Rab11;
ConcOut.R1R1_165a_int = ConcOut.R1R1_165a_Rab4+ConcOut.R1R1_165a_Rab11;

ConcOut.R1R1_plgf1_int = ConcOut.R1R1_plgf1_Rab4+ConcOut.R1R1_plgf1_Rab11;
ConcOut.R1R1_plgf2_int = ConcOut.R1R1_plgf2_Rab4+ConcOut.R1R1_plgf2_Rab11;


%% Surface percentage for Receptors

ConcOut.R1_percSurf = ConcOut.R1_surf./ConcOut.R1_total;
ConcOut.R2_percSurf = ConcOut.R2_surf./ConcOut.R2_total;
ConcOut.N1_percSurf = ConcOut.N1_surf./ConcOut.N1_total;

% remove NaN (typically due to rec levels being zero)
ConcOut.R1_percSurf(isnan(ConcOut.R1_percSurf)) = 0;
ConcOut.R2_percSurf(isnan(ConcOut.R2_percSurf)) = 0;
ConcOut.N1_percSurf(isnan(ConcOut.N1_percSurf)) = 0;
   

%% Receptors - dimeric fractions
% Parsing out the dimerization states of receptors

ConcOut.R1_dimers = observables_out(:,10)' + observables_out(:,15)' + observables_out(:,20)'; % surf, rab4, rab11
ConcOut.R1_monomers = ConcOut.R1_total - ConcOut.R1_dimers;

ConcOut.R1_dimericfraction = ConcOut.R1_dimers./ConcOut.R1_total;
ConcOut.R1_dimericfraction_surf = observables_out(:,10)'./ConcOut.R1_surf;
ConcOut.R1_dimericfraction_Rab4 = observables_out(:,15)'./ConcOut.R1_Rab4;
ConcOut.R1_dimericfraction_Rab11 = observables_out(:,20)'./ConcOut.R1_Rab11;

% remove NaN (typically due to rec levels being zero)
ConcOut.R1_dimericfraction(isnan(ConcOut.R1_dimericfraction)) = 0;
ConcOut.R1_dimericfraction_surf(isnan(ConcOut.R1_dimericfraction_surf)) = 0;
ConcOut.R1_dimericfraction_Rab4(isnan(ConcOut.R1_dimericfraction_Rab4)) = 0;
ConcOut.R1_dimericfraction_Rab11(isnan(ConcOut.R1_dimericfraction_Rab11)) = 0;

ConcOut.R2_dimers = observables_out(:,11)' + observables_out(:,16)' + observables_out(:,21)';
ConcOut.R2_monomers = ConcOut.R2_total - ConcOut.R2_dimers;

ConcOut.R2_dimericfraction = ConcOut.R2_dimers./ConcOut.R2_total;
ConcOut.R2_dimericfraction_surf = observables_out(:,11)'./ConcOut.R2_surf;
ConcOut.R2_dimericfraction_Rab4 = observables_out(:,16)'./ConcOut.R2_Rab4;
ConcOut.R2_dimericfraction_Rab11 = observables_out(:,21)'./ConcOut.R2_Rab11;

% remove NaN (typically due to rec levels being zero)
ConcOut.R2_dimericfraction(isnan(ConcOut.R2_dimericfraction)) = 0;
ConcOut.R2_dimericfraction_surf(isnan(ConcOut.R2_dimericfraction_surf)) = 0;
ConcOut.R2_dimericfraction_Rab4(isnan(ConcOut.R2_dimericfraction_Rab4)) = 0;
ConcOut.R2_dimericfraction_Rab11(isnan(ConcOut.R2_dimericfraction_Rab11)) = 0;

ConcOut.N1_dimers_surf = observables_out(:,148)' ; % added ligand-N1N1
ConcOut.N1_dimers_Rab4 = observables_out(:,149)' ; % added ligand-N1N1
ConcOut.N1_dimers_Rab11= observables_out(:,150)' ;% added ligand-N1N1

ConcOut.N1_dimers = ConcOut.N1_dimers_surf + ConcOut.N1_dimers_Rab4 + ConcOut.N1_dimers_Rab11;                      ;
ConcOut.N1_monomers = ConcOut.N1_total - ConcOut.N1_dimers;

ConcOut.N1_dimericfraction = ConcOut.N1_dimers./ConcOut.N1_total;
ConcOut.N1_dimericfraction_surf = ConcOut.N1_dimers_surf./ConcOut.N1_surf;
ConcOut.N1_dimericfraction_Rab4 = ConcOut.N1_dimers_Rab4./ConcOut.N1_Rab4;
ConcOut.N1_dimericfraction_Rab11 = ConcOut.N1_dimers_Rab11./ConcOut.N1_Rab11;

% remove NaN (typically due to rec levels being zero)
ConcOut.N1_dimericfraction(isnan(ConcOut.N1_dimericfraction)) = 0;
ConcOut.N1_dimericfraction_surf(isnan(ConcOut.N1_dimericfraction_surf)) = 0;
ConcOut.N1_dimericfraction_Rab4(isnan(ConcOut.N1_dimericfraction_Rab4)) = 0;
ConcOut.N1_dimericfraction_Rab11(isnan(ConcOut.N1_dimericfraction_Rab11)) = 0;


%% ACTIVE COMPLEXES

ConcOut.R2R2_121a_surf_ACT = observables_out(:,91)';
ConcOut.R2R2_165a_surf_ACT = observables_out(:,92)';

ConcOut.R1R1_121a_surf_ACT = observables_out(:,93)';
ConcOut.R1R1_165a_surf_ACT = observables_out(:,94)';

ConcOut.R1R1_plgf1_surf_ACT = observables_out(:,95)';
ConcOut.R1R1_plgf2_surf_ACT = observables_out(:,96)';


ConcOut.R2R2_121a_Rab4_ACT = observables_out(:,103)';
ConcOut.R2R2_165a_Rab4_ACT = observables_out(:,104)';

ConcOut.R1R1_121a_Rab4_ACT = observables_out(:,105)';
ConcOut.R1R1_165a_Rab4_ACT = observables_out(:,106)';

ConcOut.R1R1_plgf1_Rab4_ACT = observables_out(:,107)';
ConcOut.R1R1_plgf2_Rab4_ACT = observables_out(:,108)';


ConcOut.R2R2_121a_Rab11_ACT = observables_out(:,115)';
ConcOut.R2R2_165a_Rab11_ACT = observables_out(:,116)';

ConcOut.R1R1_121a_Rab11_ACT = observables_out(:,117)';
ConcOut.R1R1_165a_Rab11_ACT = observables_out(:,118)';

ConcOut.R1R1_plgf1_Rab11_ACT = observables_out(:,119)';
ConcOut.R1R1_plgf2_Rab11_ACT = observables_out(:,120)';


ConcOut.R2R2_121a_total_ACT = ConcOut.R2R2_121a_surf_ACT+ConcOut.R2R2_121a_Rab4_ACT+ConcOut.R2R2_121a_Rab11_ACT;
ConcOut.R2R2_165a_total_ACT = ConcOut.R2R2_165a_surf_ACT+ConcOut.R2R2_165a_Rab4_ACT+ConcOut.R2R2_165a_Rab11_ACT;

ConcOut.R1R1_121a_total_ACT = ConcOut.R1R1_121a_surf_ACT+ConcOut.R1R1_121a_Rab4_ACT+ConcOut.R1R1_121a_Rab11_ACT;
ConcOut.R1R1_165a_total_ACT = ConcOut.R1R1_165a_surf_ACT+ConcOut.R1R1_165a_Rab4_ACT+ConcOut.R1R1_165a_Rab11_ACT;

ConcOut.R1R1_plgf1_total_ACT = ConcOut.R1R1_plgf1_surf_ACT+ConcOut.R1R1_plgf1_Rab4_ACT+ConcOut.R1R1_plgf1_Rab11_ACT;
ConcOut.R1R1_plgf2_total_ACT = ConcOut.R1R1_plgf2_surf_ACT+ConcOut.R1R1_plgf2_Rab4_ACT+ConcOut.R1R1_plgf2_Rab11_ACT;

ConcOut.R2R2_121a_int_ACT = ConcOut.R2R2_121a_Rab4_ACT+ConcOut.R2R2_121a_Rab11_ACT;
ConcOut.R2R2_165a_int_ACT = ConcOut.R2R2_165a_Rab4_ACT+ConcOut.R2R2_165a_Rab11_ACT;

ConcOut.R1R1_121a_int_ACT = ConcOut.R1R1_121a_Rab4_ACT+ConcOut.R1R1_121a_Rab11_ACT;
ConcOut.R1R1_165a_int_ACT = ConcOut.R1R1_165a_Rab4_ACT+ConcOut.R1R1_165a_Rab11_ACT;

ConcOut.R1R1_plgf1_int_ACT = ConcOut.R1R1_plgf1_Rab4_ACT+ConcOut.R1R1_plgf1_Rab11_ACT;
ConcOut.R1R1_plgf2_int_ACT = ConcOut.R1R1_plgf2_Rab4_ACT+ConcOut.R1R1_plgf2_Rab11_ACT;


%% OVERALL ACTIVE COMPLEXES

ConcOut.R2_surf_ACT = ConcOut.R2R2_121a_surf_ACT + ConcOut.R2R2_165a_surf_ACT;
ConcOut.R1_surf_ACT = ConcOut.R1R1_121a_surf_ACT + ConcOut.R1R1_165a_surf_ACT + ConcOut.R1R1_plgf1_surf_ACT + ConcOut.R1R1_plgf2_surf_ACT;

ConcOut.R2_Rab4_ACT = ConcOut.R2R2_121a_Rab4_ACT + ConcOut.R2R2_165a_Rab4_ACT;
ConcOut.R1_Rab4_ACT = ConcOut.R1R1_121a_Rab4_ACT + ConcOut.R1R1_165a_Rab4_ACT + ConcOut.R1R1_plgf1_Rab4_ACT + ConcOut.R1R1_plgf2_Rab4_ACT;

ConcOut.R2_Rab11_ACT = ConcOut.R2R2_121a_Rab11_ACT + ConcOut.R2R2_165a_Rab11_ACT;
ConcOut.R1_Rab11_ACT = ConcOut.R1R1_121a_Rab11_ACT + ConcOut.R1R1_165a_Rab11_ACT + ConcOut.R1R1_plgf1_Rab11_ACT + ConcOut.R1R1_plgf2_Rab11_ACT;

ConcOut.R2_int_ACT = ConcOut.R2_Rab4_ACT + ConcOut.R2_Rab11_ACT;
ConcOut.R1_int_ACT = ConcOut.R1_Rab4_ACT + ConcOut.R1_Rab11_ACT;

ConcOut.R2_total_ACT = ConcOut.R2_surf_ACT + ConcOut.R2_int_ACT;
ConcOut.R1_total_ACT = ConcOut.R1_surf_ACT + ConcOut.R1_int_ACT;


%% EMPTY (unligated) RECEPTORS (could be in complex with other receptors)

ConcOut.R1_empty_surf = observables_out(:,127)';
ConcOut.R1_empty_Rab4 = observables_out(:,130)';
ConcOut.R1_empty_Rab11= observables_out(:,133)';

ConcOut.R2_empty_surf = observables_out(:,128)';
ConcOut.R2_empty_Rab4 = observables_out(:,131)';
ConcOut.R2_empty_Rab11= observables_out(:,134)';

ConcOut.N1_empty_surf = observables_out(:,129)';
ConcOut.N1_empty_Rab4 = observables_out(:,132)';
ConcOut.N1_empty_Rab11= observables_out(:,135)';



%% VEGFR1-NRP1 complexes - all
% Here we want to keep track of how many R1 and N1 are in R1-N1 complexes
% To do this, we need to separate out the true numbers of 
% the R1N1 (A), R1R1N1 (B), and N1R1R1N1 (C) complexes
% On the surface
R1N1 = observables_out(:,12)';  R1R1N1 = observables_out(:,13)'; N1R1R1N1 = observables_out(:,14)';
A = R1N1 - R1R1N1;  B = (R1R1N1 - N1R1R1N1)/2;  C = N1R1R1N1/4;
ConcOut.R1_in_R1N1_surf = A + 2*B + 2*C;
ConcOut.N1_in_R1N1_surf = A +   B + 2*C;
% In Rab4/5
R1N1 = observables_out(:,17)'; R1R1N1 = observables_out(:,18)'; N1R1R1N1 = observables_out(:,19)';
A = R1N1 - R1R1N1;  B = (R1R1N1 - N1R1R1N1)/2;  C = N1R1R1N1/4;
ConcOut.R1_in_R1N1_Rab4 = A + 2*B + 2*C;
ConcOut.N1_in_R1N1_Rab4 = A +   B + 2*C;
% In Rab11
R1N1 = observables_out(:,22)'; R1R1N1 = observables_out(:,23)'; N1R1R1N1 = observables_out(:,24)';
A = R1N1 - R1R1N1;  B = (R1R1N1 - N1R1R1N1)/2;  C = N1R1R1N1/4;
ConcOut.R1_in_R1N1_Rab11 = A + 2*B + 2*C;
ConcOut.N1_in_R1N1_Rab11 = A +   B + 2*C;
% Internal
ConcOut.R1_in_R1N1_int = ConcOut.R1_in_R1N1_Rab4 + ConcOut.R1_in_R1N1_Rab11;
ConcOut.N1_in_R1N1_int = ConcOut.N1_in_R1N1_Rab4 + ConcOut.N1_in_R1N1_Rab11;
% Total
ConcOut.R1_in_R1N1_total = ConcOut.R1_in_R1N1_surf+ConcOut.R1_in_R1N1_int;
ConcOut.N1_in_R1N1_total = ConcOut.N1_in_R1N1_surf+ConcOut.N1_in_R1N1_int;

ConcOut.R1_fractioninR1N1 = ConcOut.R1_in_R1N1_total./ConcOut.R1_total;
ConcOut.R1_fractioninR1N1_surf = ConcOut.R1_in_R1N1_surf./ConcOut.R1_surf;
ConcOut.R1_fractioninR1N1_Rab4 = ConcOut.R1_in_R1N1_Rab4./ConcOut.R1_Rab4;
ConcOut.R1_fractioninR1N1_Rab11 = ConcOut.R1_in_R1N1_Rab11./ConcOut.R1_Rab11;

ConcOut.N1_fractioninR1N1 = ConcOut.N1_in_R1N1_total./ConcOut.N1_total;
ConcOut.N1_fractioninR1N1_surf = ConcOut.N1_in_R1N1_surf./ConcOut.N1_surf;
ConcOut.N1_fractioninR1N1_Rab4 = ConcOut.N1_in_R1N1_Rab4./ConcOut.N1_Rab4;
ConcOut.N1_fractioninR1N1_Rab11 = ConcOut.N1_in_R1N1_Rab11./ConcOut.N1_Rab11;

% remove NaN (typically due to rec levels being zero)
ConcOut.R1_fractioninR1N1(isnan(ConcOut.R1_fractioninR1N1)) = 0;
ConcOut.R1_fractioninR1N1_surf(isnan(ConcOut.R1_fractioninR1N1_surf)) = 0;
ConcOut.R1_fractioninR1N1_Rab4(isnan(ConcOut.R1_fractioninR1N1_Rab4)) = 0;
ConcOut.R1_fractioninR1N1_Rab11(isnan(ConcOut.R1_fractioninR1N1_Rab11)) = 0;

ConcOut.N1_fractioninR1N1(isnan(ConcOut.N1_fractioninR1N1)) = 0;
ConcOut.N1_fractioninR1N1_surf(isnan(ConcOut.N1_fractioninR1N1_surf)) = 0;
ConcOut.N1_fractioninR1N1_Rab4(isnan(ConcOut.N1_fractioninR1N1_Rab4)) = 0;
ConcOut.N1_fractioninR1N1_Rab11(isnan(ConcOut.N1_fractioninR1N1_Rab11)) = 0;


%% VEGFR1-NRP1 complexes - active R1
% Here we want to keep track of how many R1 and N1 are in R1-N1 complexes
%  but specifically those that have active R1-R1 (i.e. with a ligand)
% To do this, we need to account for the R1R1N1 and N1R1R1N1 complexes
% On the surface
R1R1N1act   = observables_out(:,99)'  + observables_out(:,101)';  
N1R1R1N1act = observables_out(:,100)' + observables_out(:,102)'; 
ConcOut.R1_in_R1N1_surf_ACT = R1R1N1act - N1R1R1N1act/2;
ConcOut.N1_in_R1N1_surf_ACT = R1R1N1act/2;
% In Rab4/5
R1R1N1act   = observables_out(:,111)' + observables_out(:,113)';  
N1R1R1N1act = observables_out(:,112)' + observables_out(:,114)'; 
ConcOut.R1_in_R1N1_Rab4_ACT = R1R1N1act - N1R1R1N1act/2;
ConcOut.N1_in_R1N1_Rab4_ACT = R1R1N1act/2;
% In Rab11
R1R1N1act   = observables_out(:,123)' + observables_out(:,125)';  
N1R1R1N1act = observables_out(:,124)' + observables_out(:,126)'; 
ConcOut.R1_in_R1N1_Rab11_ACT = R1R1N1act - N1R1R1N1act/2;
ConcOut.N1_in_R1N1_Rab11_ACT = R1R1N1act/2;
% Internal
ConcOut.R1_in_R1N1_int_ACT = ConcOut.R1_in_R1N1_Rab4_ACT + ConcOut.R1_in_R1N1_Rab11_ACT;
ConcOut.N1_in_R1N1_int_ACT = ConcOut.N1_in_R1N1_Rab4_ACT + ConcOut.N1_in_R1N1_Rab11_ACT;
% Total
ConcOut.R1_in_R1N1_total_ACT = ConcOut.R1_in_R1N1_surf_ACT+ConcOut.R1_in_R1N1_int_ACT;
ConcOut.N1_in_R1N1_total_ACT = ConcOut.N1_in_R1N1_surf_ACT+ConcOut.N1_in_R1N1_int_ACT;

ConcOut.R1_fractioninR1N1_ACT = ConcOut.R1_in_R1N1_total_ACT./ConcOut.R1_total_ACT;
ConcOut.R1_fractioninR1N1_surf_ACT = ConcOut.R1_in_R1N1_surf_ACT./ConcOut.R1_surf_ACT;
ConcOut.R1_fractioninR1N1_Rab4_ACT = ConcOut.R1_in_R1N1_Rab4_ACT./ConcOut.R1_Rab4_ACT;
ConcOut.R1_fractioninR1N1_Rab11_ACT = ConcOut.R1_in_R1N1_Rab11_ACT./ConcOut.R1_Rab11_ACT;

% ConcOut.N1_fractioninR1N1_ACT = ConcOut.N1_in_R1N1_total_ACT./ConcOut.N1_total_ACT;
% ConcOut.N1_fractioninR1N1_surf_ACT = ConcOut.N1_in_R1N1_surf_ACT./ConcOut.N1_surf_ACT;
% ConcOut.N1_fractioninR1N1_Rab4_ACT = ConcOut.N1_in_R1N1_Rab4_ACT./ConcOut.N1_Rab4_ACT;
% ConcOut.N1_fractioninR1N1_Rab11_ACT = ConcOut.N1_in_R1N1_Rab11_ACT./ConcOut.N1_Rab11_ACT;

% remove NaN (typically due to rec levels being zero)
ConcOut.R1_fractioninR1N1_ACT(isnan(ConcOut.R1_fractioninR1N1_ACT)) = 0;
ConcOut.R1_fractioninR1N1_surf_ACT(isnan(ConcOut.R1_fractioninR1N1_surf_ACT)) = 0;
ConcOut.R1_fractioninR1N1_Rab4_ACT(isnan(ConcOut.R1_fractioninR1N1_Rab4_ACT)) = 0;
ConcOut.R1_fractioninR1N1_Rab11_ACT(isnan(ConcOut.R1_fractioninR1N1_Rab11_ACT)) = 0;

% ConcOut.N1_fractioninR1N1_ACT(isnan(ConcOut.N1_fractioninR1N1_ACT)) = 0;
% ConcOut.N1_fractioninR1N1_surf_ACT(isnan(ConcOut.N1_fractioninR1N1_surf_ACT)) = 0;
% ConcOut.N1_fractioninR1N1_Rab4_ACT(isnan(ConcOut.N1_fractioninR1N1_Rab4_ACT)) = 0;
% ConcOut.N1_fractioninR1N1_Rab11_ACT(isnan(ConcOut.N1_fractioninR1N1_Rab11_ACT)) = 0;



%% VEGFR2-NRP1 complexes - all
% Here we want to keep track of how many R2 and N1 are in R2-N1 complexes
% Which can only happen when bridged by ligands (VEGF165)
% To do this, we need to separate out the true numbers of 
% the R1N1 (A), R1R1N1 (B), and N1R1R1N1 (C) complexes
% On the surface
R2N1 = observables_out(:,151)';  R2R2N1 = observables_out(:,31)'; 
N1R2R2N1 = observables_out(:,32)'; N1R2N1 = observables_out(:,154)';
ConcOut.R2_in_R2N1_surf = R2N1 - N1R2N1/2;
ConcOut.N1_in_R2N1_surf = R2N1 - R2R2N1/2;
% In Rab4/5
R2N1 = observables_out(:,152)';  R2R2N1 = observables_out(:,43)'; 
N1R2R2N1 = observables_out(:,44)'; N1R2N1 = observables_out(:,155)';
ConcOut.R2_in_R2N1_Rab4 = R2N1 - N1R2N1/2;
ConcOut.N1_in_R2N1_Rab4 = R2N1 - R2R2N1/2;
% In Rab11
R2N1 = observables_out(:,153)';  R2R2N1 = observables_out(:,55)'; 
N1R2R2N1 = observables_out(:,56)'; N1R2N1 = observables_out(:,156)';
ConcOut.R2_in_R2N1_Rab11 = R2N1 - N1R2N1/2;
ConcOut.N1_in_R2N1_Rab11 = R2N1 - R2R2N1/2;
% Internal
ConcOut.R2_in_R2N1_int = ConcOut.R2_in_R2N1_Rab4 + ConcOut.R2_in_R2N1_Rab11;
ConcOut.N1_in_R2N1_int = ConcOut.N1_in_R2N1_Rab4 + ConcOut.N1_in_R2N1_Rab11;
% Total
ConcOut.R2_in_R2N1_total = ConcOut.R2_in_R2N1_surf+ConcOut.R2_in_R2N1_int;
ConcOut.N1_in_R2N1_total = ConcOut.N1_in_R2N1_surf+ConcOut.N1_in_R2N1_int;

ConcOut.R2_fractioninR2N1 = ConcOut.R2_in_R2N1_total./ConcOut.R2_total;
ConcOut.R2_fractioninR2N1_surf = ConcOut.R2_in_R2N1_surf./ConcOut.R2_surf;
ConcOut.R2_fractioninR2N1_Rab4 = ConcOut.R2_in_R2N1_Rab4./ConcOut.R2_Rab4;
ConcOut.R2_fractioninR2N1_Rab11 = ConcOut.R2_in_R2N1_Rab11./ConcOut.R2_Rab11;

ConcOut.N1_fractioninR2N1 = ConcOut.N1_in_R2N1_total./ConcOut.N1_total;
ConcOut.N1_fractioninR2N1_surf = ConcOut.N1_in_R2N1_surf./ConcOut.N1_surf;
ConcOut.N1_fractioninR2N1_Rab4 = ConcOut.N1_in_R2N1_Rab4./ConcOut.N1_Rab4;
ConcOut.N1_fractioninR2N1_Rab11 = ConcOut.N1_in_R2N1_Rab11./ConcOut.N1_Rab11;

% remove NaN (typically due to rec levels being zero)
ConcOut.R2_fractioninR2N1(isnan(ConcOut.R2_fractioninR2N1)) = 0;
ConcOut.R2_fractioninR2N1_surf(isnan(ConcOut.R2_fractioninR2N1_surf)) = 0;
ConcOut.R2_fractioninR2N1_Rab4(isnan(ConcOut.R2_fractioninR2N1_Rab4)) = 0;
ConcOut.R2_fractioninR2N1_Rab11(isnan(ConcOut.R2_fractioninR2N1_Rab11)) = 0;

ConcOut.N1_fractioninR2N1(isnan(ConcOut.N1_fractioninR2N1)) = 0;
ConcOut.N1_fractioninR2N1_surf(isnan(ConcOut.N1_fractioninR2N1_surf)) = 0;
ConcOut.N1_fractioninR2N1_Rab4(isnan(ConcOut.N1_fractioninR2N1_Rab4)) = 0;
ConcOut.N1_fractioninR2N1_Rab11(isnan(ConcOut.N1_fractioninR2N1_Rab11)) = 0;


%% VEGFR2-NRP1 complexes - active R2
% Here we want to keep track of how many R2 and N1 are in R2-N1 complexes
%  but specifically those that have active R2-R2 (i.e. with a ligand)
% To do this, we need to account for the R2R2N1 and N1R2R2N1 complexes
% On the surface
R2R2N1act   = observables_out(:,97)';  
N1R2R2N1act = observables_out(:,98)'; 
ConcOut.R2_in_R2N1_surf_ACT = R2R2N1act - N1R2R2N1act/2;
ConcOut.N1_in_R2N1_surf_ACT = R2R2N1act/2;
% In Rab4/5
R2R2N1act   = observables_out(:,109)';  
N1R2R2N1act = observables_out(:,110)'; 
ConcOut.R2_in_R2N1_Rab4_ACT = R2R2N1act - N1R2R2N1act/2;
ConcOut.N1_in_R2N1_Rab4_ACT = R2R2N1act/2;
% In Rab11
R2R2N1act   = observables_out(:,121)';  
N1R2R2N1act = observables_out(:,122)'; 
ConcOut.R2_in_R2N1_Rab11_ACT = R2R2N1act - N1R2R2N1act/2;
ConcOut.N1_in_R2N1_Rab11_ACT = R2R2N1act/2;
% Internal
ConcOut.R2_in_R2N1_int_ACT = ConcOut.R2_in_R2N1_Rab4_ACT + ConcOut.R2_in_R2N1_Rab11_ACT;
ConcOut.N1_in_R2N1_int_ACT = ConcOut.N1_in_R2N1_Rab4_ACT + ConcOut.N1_in_R2N1_Rab11_ACT;
% Total
ConcOut.R2_in_R2N1_total_ACT = ConcOut.R2_in_R2N1_surf_ACT+ConcOut.R2_in_R2N1_int_ACT;
ConcOut.N1_in_R2N1_total_ACT = ConcOut.N1_in_R2N1_surf_ACT+ConcOut.N1_in_R2N1_int_ACT;

ConcOut.R2_fractioninR2N1_ACT = ConcOut.R2_in_R2N1_total_ACT./ConcOut.R2_total_ACT;
ConcOut.R2_fractioninR2N1_surf_ACT = ConcOut.R2_in_R2N1_surf_ACT./ConcOut.R2_surf_ACT;
ConcOut.R2_fractioninR2N1_Rab4_ACT = ConcOut.R2_in_R2N1_Rab4_ACT./ConcOut.R2_Rab4_ACT;
ConcOut.R2_fractioninR2N1_Rab11_ACT = ConcOut.R2_in_R2N1_Rab11_ACT./ConcOut.R2_Rab11_ACT;

% ConcOut.N1_fractioninR2N1_ACT = ConcOut.N1_in_R2N1_total_ACT./ConcOut.N1_total_ACT;
% ConcOut.N1_fractioninR2N1_surf_ACT = ConcOut.N1_in_R2N1_surf_ACT./ConcOut.N1_surf_ACT;
% ConcOut.N1_fractioninR2N1_Rab4_ACT = ConcOut.N1_in_R2N1_Rab4_ACT./ConcOut.N1_Rab4_ACT;
% ConcOut.N1_fractioninR2N1_Rab11_ACT = ConcOut.N1_in_R2N1_Rab11_ACT./ConcOut.N1_Rab11_ACT;

% remove NaN (typically due to rec levels being zero)
ConcOut.R2_fractioninR2N1_ACT(isnan(ConcOut.R2_fractioninR2N1_ACT)) = 0;
ConcOut.R2_fractioninR2N1_surf_ACT(isnan(ConcOut.R2_fractioninR2N1_surf_ACT)) = 0;
ConcOut.R2_fractioninR2N1_Rab4_ACT(isnan(ConcOut.R2_fractioninR2N1_Rab4_ACT)) = 0;
ConcOut.R2_fractioninR2N1_Rab11_ACT(isnan(ConcOut.R2_fractioninR2N1_Rab11_ACT)) = 0;

% ConcOut.N1_fractioninR2N1_ACT(isnan(ConcOut.N1_fractioninR2N1_ACT)) = 0;
% ConcOut.N1_fractioninR2N1_surf_ACT(isnan(ConcOut.N1_fractioninR2N1_surf_ACT)) = 0;
% ConcOut.N1_fractioninR2N1_Rab4_ACT(isnan(ConcOut.N1_fractioninR2N1_Rab4_ACT)) = 0;
% ConcOut.N1_fractioninR2N1_Rab11_ACT(isnan(ConcOut.N1_fractioninR2N1_Rab11_ACT)) = 0;


%% LIGANDS (free ligands)

ConcOut.V165_surf = observables_out(:,61)'./4;
ConcOut.V165_Rab4 = observables_out(:,65)'./4;
ConcOut.V165_Rab11= observables_out(:,69)'./4;
ConcOut.V121_surf = observables_out(:,62)'./2;
ConcOut.V121_Rab4 = observables_out(:,66)'./2;
ConcOut.V121_Rab11= observables_out(:,70)'./2;
ConcOut.PGF1_surf = observables_out(:,63)'./2;
ConcOut.PGF1_Rab4 = observables_out(:,67)'./2;
ConcOut.PGF1_Rab11= observables_out(:,71)'./2;
ConcOut.PGF2_surf = observables_out(:,64)'./4;
ConcOut.PGF2_Rab4 = observables_out(:,68)'./4;
ConcOut.PGF2_Rab11= observables_out(:,72)'./4;

ConcOut.V165_int = ConcOut.V165_Rab4 + ConcOut.V165_Rab11 ;
ConcOut.V121_int = ConcOut.V121_Rab4 + ConcOut.V121_Rab11 ;
ConcOut.PGF1_int = ConcOut.PGF1_Rab4 + ConcOut.PGF1_Rab11 ;
ConcOut.PGF2_int = ConcOut.PGF2_Rab4 + ConcOut.PGF2_Rab11 ;

ConcOut.ligand_surf = ConcOut.V165_surf + ConcOut.V121_surf + ConcOut.PGF1_surf + ConcOut.PGF2_surf;
ConcOut.ligand_int  = ConcOut.V165_int  + ConcOut.V121_int  + ConcOut.PGF1_int  + ConcOut.PGF2_int ;
ConcOut.ligand_intperc  = ConcOut.ligand_int./ConcOut.ligand_surf*100;

%% LIGANDS (Total ligands)

ConcOut.TotV165_surf = observables_out(:,136)';
ConcOut.TotV165_Rab4 = observables_out(:,140)';
ConcOut.TotV165_Rab11= observables_out(:,144)';
ConcOut.TotV121_surf = observables_out(:,137)';
ConcOut.TotV121_Rab4 = observables_out(:,141)';
ConcOut.TotV121_Rab11= observables_out(:,145)';
ConcOut.TotPGF1_surf = observables_out(:,138)';
ConcOut.TotPGF1_Rab4 = observables_out(:,142)';
ConcOut.TotPGF1_Rab11= observables_out(:,146)';
ConcOut.TotPGF2_surf = observables_out(:,139)';
ConcOut.TotPGF2_Rab4 = observables_out(:,143)';
ConcOut.TotPGF2_Rab11= observables_out(:,147)';

ConcOut.TotV165_int = ConcOut.TotV165_Rab4 + ConcOut.TotV165_Rab11 ;
ConcOut.TotV121_int = ConcOut.TotV121_Rab4 + ConcOut.TotV121_Rab11 ;
ConcOut.TotPGF1_int = ConcOut.TotPGF1_Rab4 + ConcOut.TotPGF1_Rab11 ;
ConcOut.TotPGF2_int = ConcOut.TotPGF2_Rab4 + ConcOut.TotPGF2_Rab11 ;

%% LIGANDS (Bound ligands)

ConcOut.BndV165_surf = ConcOut.TotV165_surf - ConcOut.V165_surf;
ConcOut.BndV165_Rab4 = ConcOut.TotV165_Rab4 - ConcOut.V165_Rab4;
ConcOut.BndV165_Rab11= ConcOut.TotV165_Rab11 - ConcOut.V165_Rab11;
ConcOut.BndV121_surf = ConcOut.TotV121_surf - ConcOut.V121_surf;
ConcOut.BndV121_Rab4 = ConcOut.TotV121_Rab4 - ConcOut.V121_Rab4;
ConcOut.BndV121_Rab11= ConcOut.TotV121_Rab11 - ConcOut.V121_Rab11;
ConcOut.BndPGF1_surf = ConcOut.TotPGF1_surf - ConcOut.PGF1_surf;
ConcOut.BndPGF1_Rab4 = ConcOut.TotPGF1_Rab4 - ConcOut.PGF1_Rab4;
ConcOut.BndPGF1_Rab11= ConcOut.TotPGF1_Rab11 - ConcOut.PGF1_Rab11;
ConcOut.BndPGF2_surf = ConcOut.TotPGF2_surf - ConcOut.PGF2_surf;
ConcOut.BndPGF2_Rab4 = ConcOut.TotPGF2_Rab4 - ConcOut.PGF2_Rab4;
ConcOut.BndPGF2_Rab11= ConcOut.TotPGF2_Rab11 - ConcOut.PGF2_Rab11;

ConcOut.BndV165_int = ConcOut.BndV165_Rab4 + ConcOut.BndV165_Rab11 ;
ConcOut.BndV121_int = ConcOut.BndV121_Rab4 + ConcOut.BndV121_Rab11 ;
ConcOut.BndPGF1_int = ConcOut.BndPGF1_Rab4 + ConcOut.BndPGF1_Rab11 ;
ConcOut.BndPGF2_int = ConcOut.BndPGF2_Rab4 + ConcOut.BndPGF2_Rab11 ;

ConcOut.BndLigand_surf = ConcOut.BndV165_surf + ConcOut.BndV121_surf + ConcOut.BndPGF1_surf + ConcOut.BndPGF2_surf;
ConcOut.BndLigand_int  = ConcOut.BndV165_int  + ConcOut.BndV121_int  + ConcOut.BndPGF1_int  + ConcOut.BndPGF2_int ;

ConcOut.BndInactiveLig_surf = ConcOut.BndLigand_surf - ConcOut.R1_surf_ACT/2.0 - ConcOut.R2_surf_ACT/2.0;
ConcOut.BndInactiveLig_int  = ConcOut.BndLigand_int  - ConcOut.R1_int_ACT/2.0  - ConcOut.R2_int_ACT/2.0;

ConcOut.BndInactiveLig_surfperc = ConcOut.BndInactiveLig_surf./ConcOut.BndLigand_surf*100;
ConcOut.BndInactiveLig_intperc  = ConcOut.BndInactiveLig_int./ConcOut.BndLigand_int*100;






