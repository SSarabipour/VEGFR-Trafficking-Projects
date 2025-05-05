function ConcOut = CalcOutputs(observables_out, suppressligands)
% Identifying the key observable metrics, and assembling them into a structure
%
% INPUTS: observables_out = a set of observable metrics, as output by the simulation code
%         
% OUTPUTS: ConcOut = stucture of key metrics - Total Receptors 
%

%% Receptors - whole cell and localized
% In all cases what's being counted are total monomers, including monomers
% involved in dimers and higher order complexes. 

ConcOut.R1_surf  = observables_out(:,4)';
ConcOut.R2_surf  = observables_out(:,5)';
ConcOut.N1_surf  = observables_out(:,6)';

ConcOut.R1_Rab4  = observables_out(:,7)';
ConcOut.R2_Rab4  = observables_out(:,8)';
ConcOut.N1_Rab4  = observables_out(:,9)';

ConcOut.R1_Rab11 = observables_out(:,10)';
ConcOut.R2_Rab11 = observables_out(:,11)';
ConcOut.N1_Rab11 = observables_out(:,12)';

ConcOut.R1_int   =                 ConcOut.R1_Rab4+ConcOut.R1_Rab11;
ConcOut.R2_int   =                 ConcOut.R2_Rab4+ConcOut.R2_Rab11;
ConcOut.N1_int   =                 ConcOut.N1_Rab4+ConcOut.N1_Rab11;

ConcOut.R1_Rab7  = observables_out(:,13)'; 
ConcOut.R2_Rab7  = observables_out(:,14)';
ConcOut.N1_Rab7  = observables_out(:,15)';

ConcOut.R1_total = ConcOut.R1_surf+ConcOut.R1_int;
ConcOut.R2_total = ConcOut.R2_surf+ConcOut.R2_int;
ConcOut.N1_total = ConcOut.N1_surf+ConcOut.N1_int;


%% Receptors - dimeric fractions
% Parsing out the dimerization states of receptors

ConcOut.R1_dimers = observables_out(:,21)' + observables_out(:,26)' + observables_out(:,31)';
ConcOut.R1_monomers = ConcOut.R1_total - ConcOut.R1_dimers;
if ConcOut.R1_total == 0
    ConcOut.R1_dimericfraction = 0;
else
    ConcOut.R1_dimericfraction = ConcOut.R1_dimers./ConcOut.R1_total;
end
if ConcOut.R1_surf == 0
    ConcOut.R1_dimericfraction_surf = 0;
else
    ConcOut.R1_dimericfraction_surf = observables_out(:,21)'./ConcOut.R1_surf;
end
if ConcOut.R1_Rab4 == 0
    ConcOut.R1_dimericfraction_Rab4 = 0;
else
    ConcOut.R1_dimericfraction_Rab4 = observables_out(:,26)'./ConcOut.R1_Rab4;
end
if ConcOut.R1_Rab11 == 0
    ConcOut.R1_dimericfraction_Rab11 = 0;
else
    ConcOut.R1_dimericfraction_Rab11 = observables_out(:,31)'./ConcOut.R1_Rab11;
end

ConcOut.R2_dimers = observables_out(:,22)' + observables_out(:,27)' + observables_out(:,32)';
ConcOut.R2_monomers = ConcOut.R2_total - ConcOut.R2_dimers;
if ConcOut.R2_total == 0
    ConcOut.R2_dimericfraction = 0;
else
    ConcOut.R2_dimericfraction = ConcOut.R2_dimers./ConcOut.R2_total;
end
if ConcOut.R2_surf == 0
    ConcOut.R2_dimericfraction_surf = 0;
else
    ConcOut.R2_dimericfraction_surf = observables_out(:,22)'./ConcOut.R2_surf;
end
if ConcOut.R2_Rab4 == 0
    ConcOut.R2_dimericfraction_Rab4 = 0;
else
    ConcOut.R2_dimericfraction_Rab4 = observables_out(:,27)'./ConcOut.R2_Rab4;
end
if ConcOut.R2_Rab11 == 0
    ConcOut.R2_dimericfraction_Rab11 = 0;
else
    ConcOut.R2_dimericfraction_Rab11 = observables_out(:,32)'./ConcOut.R2_Rab11;
end

ConcOut.N1_dimers = observables_out(:,25)' + observables_out(:,30)' + observables_out(:,35)';
ConcOut.N1_monomers = ConcOut.N1_total - ConcOut.N1_dimers;
if ConcOut.N1_total == 0
    ConcOut.N1_dimericfraction = 0;
else
    ConcOut.N1_dimericfraction = ConcOut.N1_dimers./ConcOut.N1_total;
end
if ConcOut.N1_surf == 0
    ConcOut.N1_dimericfraction_surf = 0;
else
    ConcOut.N1_dimericfraction_surf = observables_out(:,25)'./ConcOut.N1_surf;
end
if ConcOut.N1_Rab4 == 0
    ConcOut.N1_dimericfraction_Rab4 = 0;
else
    ConcOut.N1_dimericfraction_Rab4 = observables_out(:,30)'./ConcOut.N1_Rab4;
end
if ConcOut.N1_Rab11 == 0
    ConcOut.N1_dimericfraction_Rab11 = 0;
else
    ConcOut.N1_dimericfraction_Rab11 = observables_out(:,35)'./ConcOut.N1_Rab11;
end


%% VEGFR1-NRP1 complexes
% Here we want to keep track of how many R1 and N1 are in R1-N1 complexes
% To do this, we need to separate out the true numbers of 
% the R1N1 (A), R1R1N1 (B), and N1R1R1N1 (C) complexes
% On the surface
R1N1 = observables_out(:,23)'; R1R1N1 = observables_out(:,24)'; N1R1R1N1 = observables_out(:,25)';
A = R1N1 - R1R1N1;  B = (R1R1N1 - N1R1R1N1)/2;  C = N1R1R1N1/4;
ConcOut.R1_in_R1N1_surf = A + 2*B + 2*C;
ConcOut.N1_in_R1N1_surf = A +   B + 2*C;
% In Rab4/5
R1N1 = observables_out(:,28)'; R1R1N1 = observables_out(:,29)'; N1R1R1N1 = observables_out(:,30)';
A = R1N1 - R1R1N1;  B = (R1R1N1 - N1R1R1N1)/2;  C = N1R1R1N1/4;
ConcOut.R1_in_R1N1_Rab4 = A + 2*B + 2*C;
ConcOut.N1_in_R1N1_Rab4 = A +   B + 2*C;
% In Rab11
R1N1 = observables_out(:,33)'; R1R1N1 = observables_out(:,34)'; N1R1R1N1 = observables_out(:,35)';
A = R1N1 - R1R1N1;  B = (R1R1N1 - N1R1R1N1)/2;  C = N1R1R1N1/4;
ConcOut.R1_in_R1N1_Rab11 = A + 2*B + 2*C;
ConcOut.N1_in_R1N1_Rab11 = A +   B + 2*C;
% Internal
ConcOut.R1_in_R1N1_int = ConcOut.R1_in_R1N1_Rab4 + ConcOut.R1_in_R1N1_Rab11;
ConcOut.N1_in_R1N1_int = ConcOut.N1_in_R1N1_Rab4 + ConcOut.N1_in_R1N1_Rab11;
% Total
ConcOut.R1_in_R1N1_total = ConcOut.R1_in_R1N1_surf+ConcOut.R1_in_R1N1_int;
ConcOut.N1_in_R1N1_total = ConcOut.N1_in_R1N1_surf+ConcOut.N1_in_R1N1_int;

if ConcOut.R1_total == 0
    ConcOut.R1_fractioninR1N1 = 0;
else
    ConcOut.R1_fractioninR1N1 = ConcOut.R1_in_R1N1_total./ConcOut.R1_total;
end
if ConcOut.R1_surf == 0
    ConcOut.R1_fractioninR1N1_surf = 0;
else
    ConcOut.R1_fractioninR1N1_surf = ConcOut.R1_in_R1N1_surf./ConcOut.R1_surf;
end
if ConcOut.R1_Rab4 == 0
    ConcOut.R1_fractioninR1N1_Rab4 = 0;
else
    ConcOut.R1_fractioninR1N1_Rab4 = ConcOut.R1_in_R1N1_Rab4./ConcOut.R1_Rab4;
end
if ConcOut.R1_Rab11 == 0
    ConcOut.R1_fractioninR1N1_Rab11 = 0;
else
    ConcOut.R1_fractioninR1N1_Rab11 = ConcOut.R1_in_R1N1_Rab11./ConcOut.R1_Rab11;
end

if ConcOut.N1_total == 0
    ConcOut.N1_fractioninR1N1 = 0;
else
    ConcOut.N1_fractioninR1N1 = ConcOut.N1_in_R1N1_total./ConcOut.N1_total;
end
if ConcOut.N1_surf == 0
    ConcOut.N1_fractioninR1N1_surf = 0;
else
    ConcOut.N1_fractioninR1N1_surf = ConcOut.N1_in_R1N1_surf./ConcOut.N1_surf;
end
if ConcOut.N1_Rab4 == 0
    ConcOut.N1_fractioninR1N1_Rab4 = 0;
else
    ConcOut.N1_fractioninR1N1_Rab4 = ConcOut.N1_in_R1N1_Rab4./ConcOut.N1_Rab4;
end
if ConcOut.N1_Rab11 == 0
    ConcOut.N1_fractioninR1N1_Rab11 = 0;
else
    ConcOut.N1_fractioninR1N1_Rab11 = ConcOut.N1_in_R1N1_Rab11./ConcOut.N1_Rab11;
end

%% Surface percentage for Receptors

if ConcOut.R1_total == 0
    ConcOut.R1_percSurf = 0;
else
    ConcOut.R1_percSurf = ConcOut.R1_surf./ConcOut.R1_total;
end
if ConcOut.R2_total == 0
    ConcOut.R2_percSurf = 0;
else
    ConcOut.R2_percSurf = ConcOut.R2_surf./ConcOut.R2_total;
end
if ConcOut.N1_total == 0
    ConcOut.N1_percSurf = 0;
else
    ConcOut.N1_percSurf = ConcOut.N1_surf./ConcOut.N1_total;
end

end
