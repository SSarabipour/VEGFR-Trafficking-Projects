function [R1rates,R2rates,N1rates] = CalcRates(ConcOut,parameters,perturbs)
% Using the observed levels and parameter values, quantify the rates of
% transport for VEGFR1, VEGFR2, and NRP1
%
% INPUTS: ConcOut = key metrics of concentration, output from simulations
%         parameters = model parameter values
%         parameters = perturbation values (?)
%         
% OUTPUTS: Xrates= rates of transport of receptor X
%

%% Receptors - whole cell and localized
outnames = fieldnames(ConcOut);
tpts = length(ConcOut.(outnames{1}));

% names for each transport process
ratenames = {'Surf_in_Rab11','Surf_in_Rab4','Surf_in_Synth','Surf_out_Rab4',...
    'Rab4_in_Surf','Rab4_out_Rab11','Rab4_out_Surf','Rab4_out_Rab7',...
    'Rab11_in_Rab4','Rab11_out_Surf',...
    'Surf_net','Rab4_net','Rab11_net'};

for i = 1:length(ratenames)
    R1rates.(ratenames{i}) = 0;
end
R2rates = R1rates; N1rates = R1rates;


R1rates.Surf_in_Synth(1) = parameters.kR1prod;
R1rates.Surf_in_Synth(2:length(ConcOut.R1_surf)) = parameters.kR1prod * perturbs.R1reductionByCHX;
R1rates.Surf_in_Rab4  = parameters.kR1Rab4a  * (ConcOut.R1_Rab4  - ConcOut.R1_in_R1N1_Rab4)     + parameters.kR1N1Rab4a  * ConcOut.R1_in_R1N1_Rab4;
R1rates.Surf_in_Rab11 = parameters.kR1Rab11a * (ConcOut.R1_Rab11 - ConcOut.R1_in_R1N1_Rab11)    + parameters.kR1N1Rab11a  * ConcOut.R1_in_R1N1_Rab11;
R1rates.Surf_out_Rab4 = parameters.kR1Rab5a  * (ConcOut.R1_surf  - ConcOut.R1_in_R1N1_surf)  + parameters.kR1N1Rab5a  * ConcOut.R1_in_R1N1_surf;

R1rates.Rab4_in_Surf  = parameters.kR1Rab5a     * (ConcOut.R1_surf - ConcOut.R1_in_R1N1_surf)  + parameters.kR1N1Rab5a      * ConcOut.R1_in_R1N1_surf;
R1rates.Rab4_out_Rab7 = parameters.kR1Rab4at7a  * (ConcOut.R1_Rab4 - ConcOut.R1_in_R1N1_Rab4)  + parameters.kR1N1Rab4at7a   * ConcOut.R1_in_R1N1_Rab4;
R1rates.Rab4_out_Surf = parameters.kR1Rab4a     * (ConcOut.R1_Rab4 - ConcOut.R1_in_R1N1_Rab4)  + parameters.kR1N1Rab4a      * ConcOut.R1_in_R1N1_Rab4;
R1rates.Rab4_out_Rab11= parameters.kR1Rab4at11a * (ConcOut.R1_Rab4 - ConcOut.R1_in_R1N1_Rab4)  + parameters.kR1N1Rab4at11a  * ConcOut.R1_in_R1N1_Rab4;

R1rates.Rab11_in_Rab4 = parameters.kR1Rab4at11a * (ConcOut.R1_Rab4  - ConcOut.R1_in_R1N1_Rab4)   + parameters.kR1N1Rab4at11a  * ConcOut.R1_in_R1N1_Rab4;
R1rates.Rab11_out_Surf= parameters.kR1Rab11a    * (ConcOut.R1_Rab11 - ConcOut.R1_in_R1N1_Rab11)  + parameters.kR1N1Rab11a     * ConcOut.R1_in_R1N1_Rab11;

R1rates.Surf_net = R1rates.Surf_in_Synth + R1rates.Surf_in_Rab4  + R1rates.Surf_in_Rab11 - R1rates.Surf_out_Rab4;
R1rates.Rab4_net = R1rates.Rab4_in_Surf  - R1rates.Rab4_out_Rab7 - R1rates.Rab4_out_Surf - R1rates.Rab4_out_Rab11;
R1rates.Rab11_net= R1rates.Rab11_in_Rab4 - R1rates.Rab11_out_Surf;


R2rates.Surf_in_Synth(1) = parameters.kR2prod;
R2rates.Surf_in_Synth(2:length(ConcOut.R2_surf)) = parameters.kR2prod * perturbs.R2reductionByCHX;
R2rates.Surf_in_Rab4  = parameters.kR2Rab4a  * ConcOut.R2_Rab4;
R2rates.Surf_in_Rab11 = parameters.kR2Rab11a * ConcOut.R2_Rab11;
R2rates.Surf_out_Rab4 = parameters.kR2Rab5a  * ConcOut.R2_surf;

R2rates.Rab4_in_Surf  = parameters.kR2Rab5a     * ConcOut.R2_surf;
R2rates.Rab4_out_Rab7 = parameters.kR2Rab4at7a  * ConcOut.R2_Rab4;
R2rates.Rab4_out_Surf = parameters.kR2Rab4a     * ConcOut.R2_Rab4;
R2rates.Rab4_out_Rab11= parameters.kR2Rab4at11a * ConcOut.R2_Rab4;

R2rates.Rab11_in_Rab4 = parameters.kR2Rab4at11a * ConcOut.R2_Rab4;
R2rates.Rab11_out_Surf= parameters.kR2Rab11a    * ConcOut.R2_Rab11;

R2rates.Surf_net = R2rates.Surf_in_Synth + R2rates.Surf_in_Rab4  + R2rates.Surf_in_Rab11 - R2rates.Surf_out_Rab4;
R2rates.Rab4_net = R2rates.Rab4_in_Surf  - R2rates.Rab4_out_Rab7 - R2rates.Rab4_out_Surf - R2rates.Rab4_out_Rab11;
R2rates.Rab11_net= R2rates.Rab11_in_Rab4 - R2rates.Rab11_out_Surf;


N1rates.Surf_in_Synth(1) = parameters.kN1prod;
N1rates.Surf_in_Synth(2:length(ConcOut.N1_surf)) = parameters.kN1prod * perturbs.N1reductionByCHX;
N1rates.Surf_in_Rab4  = parameters.kN1Rab4a  * (ConcOut.N1_Rab4  - ConcOut.N1_in_R1N1_Rab4)  + parameters.kR1N1Rab4a  * ConcOut.N1_in_R1N1_Rab4;
N1rates.Surf_in_Rab11 = parameters.kN1Rab11a * (ConcOut.N1_Rab11 - ConcOut.N1_in_R1N1_Rab11) + parameters.kR1N1Rab11a * ConcOut.N1_in_R1N1_Rab11;
N1rates.Surf_out_Rab4 = parameters.kN1Rab5a  * (ConcOut.N1_surf  - ConcOut.N1_in_R1N1_surf)  + parameters.kR1N1Rab5a  * ConcOut.N1_in_R1N1_surf;

N1rates.Rab4_in_Surf  = parameters.kN1Rab5a     * (ConcOut.N1_surf - ConcOut.N1_in_R1N1_surf)  + parameters.kR1N1Rab5a      * ConcOut.N1_in_R1N1_surf;
N1rates.Rab4_out_Rab7 = parameters.kN1Rab4at7a  * (ConcOut.N1_Rab4 - ConcOut.N1_in_R1N1_Rab4)  + parameters.kR1N1Rab4at7a   * ConcOut.N1_in_R1N1_Rab4;
N1rates.Rab4_out_Surf = parameters.kN1Rab4a     * (ConcOut.N1_Rab4 - ConcOut.N1_in_R1N1_Rab4)  + parameters.kR1N1Rab4a      * ConcOut.N1_in_R1N1_Rab4;
N1rates.Rab4_out_Rab11= parameters.kN1Rab4at11a * (ConcOut.N1_Rab4 - ConcOut.N1_in_R1N1_Rab4)  + parameters.kR1N1Rab4at11a  * ConcOut.N1_in_R1N1_Rab4;

N1rates.Rab11_in_Rab4 = parameters.kN1Rab4at11a * (ConcOut.N1_Rab4  - ConcOut.N1_in_R1N1_Rab4)   + parameters.kR1N1Rab4at11a  * ConcOut.N1_in_R1N1_Rab4;
N1rates.Rab11_out_Surf= parameters.kN1Rab11a    * (ConcOut.N1_Rab11 - ConcOut.N1_in_R1N1_Rab11)  + parameters.kR1N1Rab11a     * ConcOut.N1_in_R1N1_Rab11;

N1rates.Surf_net = N1rates.Surf_in_Synth + N1rates.Surf_in_Rab4  + N1rates.Surf_in_Rab11 - N1rates.Surf_out_Rab4;
N1rates.Rab4_net = N1rates.Rab4_in_Surf  - N1rates.Rab4_out_Rab7 - N1rates.Rab4_out_Surf - N1rates.Rab4_out_Rab11;
N1rates.Rab11_net= N1rates.Rab11_in_Rab4 - N1rates.Rab11_out_Surf;

end
