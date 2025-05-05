function [R1rates,R2rates,N1rates] = CalcRatesL(ConcOut,parameters) % ,perturbs
% Using the observed levels and parameter values, quantify the rates of
% transport for VEGFR1, VEGFR2, and NRP1
%
% INPUTS: ConcOut = key metrics of concentration, output from simulations
%         parameters = model parameter values
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


% VEGFR1 trafficking
  
ConcOut.R1_surf_inact = ConcOut.R1_surf - ConcOut.R1_surf_ACT;
ConcOut.R1_in_R1N1_surf_inact = ConcOut.R1_in_R1N1_surf - ConcOut.R1_in_R1N1_surf_ACT;
ConcOut.R1_Rab4_inact = ConcOut.R1_Rab4 - ConcOut.R1_Rab4_ACT;
ConcOut.R1_in_R1N1_Rab4_inact = ConcOut.R1_in_R1N1_Rab4 - ConcOut.R1_in_R1N1_Rab4_ACT;
ConcOut.R1_Rab11_inact = ConcOut.R1_Rab11 - ConcOut.R1_Rab11_ACT;
ConcOut.R1_in_R1N1_Rab11_inact = ConcOut.R1_in_R1N1_Rab11 - ConcOut.R1_in_R1N1_Rab11_ACT;

internalization = parameters.kR1Rab5a  * (ConcOut.R1_surf_inact  - ConcOut.R1_in_R1N1_surf_inact)    + parameters.kR1N1Rab5a  * ConcOut.R1_in_R1N1_surf_inact ...
                + parameters.kVegfR1Rab5a  * (ConcOut.R1_surf_ACT  - ConcOut.R1_in_R1N1_surf_ACT)    + parameters.kVegfR1N1Rab5a  * ConcOut.R1_in_R1N1_surf_ACT;
degradation     = parameters.kR1Rab4at7a  * (ConcOut.R1_Rab4_inact - ConcOut.R1_in_R1N1_Rab4_inact)  + parameters.kR1N1Rab4at7a   * ConcOut.R1_in_R1N1_Rab4_inact  ...
                + parameters.kVegfR1Rab4at7a  * (ConcOut.R1_Rab4_ACT - ConcOut.R1_in_R1N1_Rab4_ACT)  + parameters.kVegfR1N1Rab4at7a   * ConcOut.R1_in_R1N1_Rab4_ACT;
recyclingRab4   = parameters.kR1Rab4a  * (ConcOut.R1_Rab4_inact  - ConcOut.R1_in_R1N1_Rab4_inact)    + parameters.kR1N1Rab4a  * ConcOut.R1_in_R1N1_Rab4_inact ...
                + parameters.kVegfR1Rab4a  * (ConcOut.R1_Rab4_ACT  - ConcOut.R1_in_R1N1_Rab4_ACT)    + parameters.kVegfR1N1Rab4a  * ConcOut.R1_in_R1N1_Rab4_ACT;
transferRab411  = parameters.kR1Rab4at11a  * (ConcOut.R1_Rab4_inact  - ConcOut.R1_in_R1N1_Rab4_inact)+ parameters.kR1N1Rab4at11a  * ConcOut.R1_in_R1N1_Rab4_inact ...
                + parameters.kVegfR1Rab4at11a  * (ConcOut.R1_Rab4_ACT  - ConcOut.R1_in_R1N1_Rab4_ACT)+ parameters.kVegfR1N1Rab4at11a  * ConcOut.R1_in_R1N1_Rab4_ACT;
recyclingRab11  = parameters.kR1Rab11a  * (ConcOut.R1_Rab11_inact  - ConcOut.R1_in_R1N1_Rab11_inact) + parameters.kR1N1Rab11a  * ConcOut.R1_in_R1N1_Rab11_inact ...
                + parameters.kVegfR1Rab11a  * (ConcOut.R1_Rab11_ACT  - ConcOut.R1_in_R1N1_Rab11_ACT) + parameters.kVegfR1N1Rab11a  * ConcOut.R1_in_R1N1_Rab11_ACT;
% NB would have to reformat these if there are differential rate constants values for
% trafficking of VEGFR1 with PlGF bound vs VEGF bound (currently there are not)

R1rates.Surf_in_Synth(1:length(ConcOut.R1_surf)) = parameters.kR1prod;
% R1rates.Surf_in_Synth(2:length(ConcOut.R1_surf)) = parameters.kR1prod * perturbs.R1reductionByCHX;
R1rates.Surf_in_Rab4  = recyclingRab4; 
R1rates.Surf_in_Rab11 = recyclingRab11; 
R1rates.Surf_out_Rab4 = internalization; 

R1rates.Rab4_in_Surf  = internalization; 
R1rates.Rab4_out_Rab7 = degradation; 
R1rates.Rab4_out_Surf = recyclingRab4; 
R1rates.Rab4_out_Rab11= transferRab411; 

R1rates.Rab11_in_Rab4 = transferRab411; 
R1rates.Rab11_out_Surf= recyclingRab11; 

R1rates.Surf_net = R1rates.Surf_in_Synth + R1rates.Surf_in_Rab4  + R1rates.Surf_in_Rab11 - R1rates.Surf_out_Rab4;
R1rates.Rab4_net = R1rates.Rab4_in_Surf  - R1rates.Rab4_out_Rab7 - R1rates.Rab4_out_Surf - R1rates.Rab4_out_Rab11;
R1rates.Rab11_net= R1rates.Rab11_in_Rab4 - R1rates.Rab11_out_Surf;


% VEGFR2 trafficking

ConcOut.R2_surf_inact = ConcOut.R2_surf - ConcOut.R2_surf_ACT;
ConcOut.R2_in_R2N1_surf_inact = ConcOut.R2_in_R2N1_surf - ConcOut.R2_in_R2N1_surf_ACT;
ConcOut.R2_Rab4_inact = ConcOut.R2_Rab4 - ConcOut.R2_Rab4_ACT;
ConcOut.R2_in_R2N1_Rab4_inact = ConcOut.R2_in_R2N1_Rab4 - ConcOut.R2_in_R2N1_Rab4_ACT;
ConcOut.R2_Rab11_inact = ConcOut.R2_Rab11 - ConcOut.R2_Rab11_ACT;
ConcOut.R2_in_R2N1_Rab11_inact = ConcOut.R2_in_R2N1_Rab11 - ConcOut.R2_in_R2N1_Rab11_ACT;

internalization = parameters.kR2Rab5a  * (ConcOut.R2_surf_inact  - ConcOut.R2_in_R2N1_surf_inact)    + parameters.kR2N1Rab5a  * ConcOut.R2_in_R2N1_surf_inact ...
                + parameters.kVegfR2Rab5a  * (ConcOut.R2_surf_ACT  - ConcOut.R2_in_R2N1_surf_ACT)    + parameters.kVegfR2N1Rab5a  * ConcOut.R2_in_R2N1_surf_ACT;
degradation     = parameters.kR2Rab4at7a  * (ConcOut.R2_Rab4_inact - ConcOut.R2_in_R2N1_Rab4_inact)  + parameters.kR2N1Rab4at7a   * ConcOut.R2_in_R2N1_Rab4_inact ...
                + parameters.kVegfR2Rab4at7a  * (ConcOut.R2_Rab4_ACT - ConcOut.R2_in_R2N1_Rab4_ACT)  + parameters.kVegfR2N1Rab4at7a   * ConcOut.R2_in_R2N1_Rab4_ACT;
recyclingRab4   = parameters.kR2Rab4a  * (ConcOut.R2_Rab4_inact  - ConcOut.R2_in_R2N1_Rab4_inact)    + parameters.kR2N1Rab4a  * ConcOut.R2_in_R2N1_Rab4_inact ...
                + parameters.kVegfR2Rab4a  * (ConcOut.R2_Rab4_ACT  - ConcOut.R2_in_R2N1_Rab4_ACT)    + parameters.kVegfR2N1Rab4a  * ConcOut.R2_in_R2N1_Rab4_ACT;
transferRab411  = parameters.kR2Rab4at11a  * (ConcOut.R2_Rab4_inact  - ConcOut.R2_in_R2N1_Rab4_inact)+ parameters.kR2N1Rab4at11a  * ConcOut.R2_in_R2N1_Rab4_inact ...
                + parameters.kVegfR2Rab4at11a  * (ConcOut.R2_Rab4_ACT  - ConcOut.R2_in_R2N1_Rab4_ACT)+ parameters.kVegfR2N1Rab4at11a  * ConcOut.R2_in_R2N1_Rab4_ACT;
recyclingRab11  = parameters.kR2Rab11a  * (ConcOut.R2_Rab11_inact  - ConcOut.R2_in_R2N1_Rab11_inact) + parameters.kR2N1Rab11a  * ConcOut.R2_in_R2N1_Rab11_inact ...
                + parameters.kVegfR2Rab11a  * (ConcOut.R2_Rab11_ACT  - ConcOut.R2_in_R2N1_Rab11_ACT) + parameters.kVegfR2N1Rab11a  * ConcOut.R2_in_R2N1_Rab11_ACT;

R2rates.Surf_in_Synth(1:length(ConcOut.R2_surf)) = parameters.kR2prod;
% R2rates.Surf_in_Synth(2:length(ConcOut.R2_surf)) = parameters.kR2prod * perturbs.R2reductionByCHX;
R2rates.Surf_in_Rab4  = recyclingRab4;
R2rates.Surf_in_Rab11 = recyclingRab11;
R2rates.Surf_out_Rab4 = internalization;

R2rates.Rab4_in_Surf  = internalization;
R2rates.Rab4_out_Rab7 = degradation;
R2rates.Rab4_out_Surf = recyclingRab4;
R2rates.Rab4_out_Rab11= transferRab411;

R2rates.Rab11_in_Rab4 = transferRab411;
R2rates.Rab11_out_Surf= recyclingRab11;

R2rates.Surf_net = R2rates.Surf_in_Synth + R2rates.Surf_in_Rab4  + R2rates.Surf_in_Rab11 - R2rates.Surf_out_Rab4;
R2rates.Rab4_net = R2rates.Rab4_in_Surf  - R2rates.Rab4_out_Rab7 - R2rates.Rab4_out_Surf - R2rates.Rab4_out_Rab11;
R2rates.Rab11_net= R2rates.Rab11_in_Rab4 - R2rates.Rab11_out_Surf;


% NRP1 trafficking

ConcOut.N1_in_R1N1_surf_inact = ConcOut.N1_in_R1N1_surf - ConcOut.N1_in_R1N1_surf_ACT;
ConcOut.N1_in_R1N1_Rab4_inact = ConcOut.N1_in_R1N1_Rab4 - ConcOut.N1_in_R1N1_Rab4_ACT;
ConcOut.N1_in_R1N1_Rab11_inact = ConcOut.N1_in_R1N1_Rab11 - ConcOut.N1_in_R1N1_Rab11_ACT;
ConcOut.N1_in_R2N1_surf_inact = ConcOut.N1_in_R2N1_surf - ConcOut.N1_in_R2N1_surf_ACT;
ConcOut.N1_in_R2N1_Rab4_inact = ConcOut.N1_in_R2N1_Rab4 - ConcOut.N1_in_R2N1_Rab4_ACT;
ConcOut.N1_in_R2N1_Rab11_inact = ConcOut.N1_in_R2N1_Rab11 - ConcOut.N1_in_R2N1_Rab11_ACT;

internalization = parameters.kN1Rab5a  * (ConcOut.N1_surf  - ConcOut.N1_in_R1N1_surf_inact - ConcOut.N1_in_R1N1_surf_ACT - ConcOut.N1_in_R2N1_surf_inact - ConcOut.N1_in_R2N1_surf_ACT) ...   
                + parameters.kR1N1Rab5a  * ConcOut.N1_in_R1N1_surf_inact + parameters.kVegfR1N1Rab5a  * ConcOut.N1_in_R1N1_surf_ACT ...
                + parameters.kR2N1Rab5a  * ConcOut.N1_in_R2N1_surf_inact + parameters.kVegfR2N1Rab5a  * ConcOut.N1_in_R2N1_surf_ACT;
degradation     = parameters.kN1Rab4at7a  * (ConcOut.N1_Rab4  - ConcOut.N1_in_R1N1_Rab4_inact - ConcOut.N1_in_R1N1_Rab4_ACT - ConcOut.N1_in_R2N1_Rab4_inact - ConcOut.N1_in_R2N1_Rab4_ACT) ...   
                + parameters.kR1N1Rab4at7a  * ConcOut.N1_in_R1N1_Rab4_inact + parameters.kVegfR1N1Rab4at7a  * ConcOut.N1_in_R1N1_Rab4_ACT ...
                + parameters.kR2N1Rab4at7a  * ConcOut.N1_in_R2N1_Rab4_inact + parameters.kVegfR2N1Rab4at7a  * ConcOut.N1_in_R2N1_Rab4_ACT;
recyclingRab4   = parameters.kN1Rab4a  * (ConcOut.N1_Rab4  - ConcOut.N1_in_R1N1_Rab4_inact - ConcOut.N1_in_R1N1_Rab4_ACT - ConcOut.N1_in_R2N1_Rab4_inact - ConcOut.N1_in_R2N1_Rab4_ACT) ...   
                + parameters.kR1N1Rab4a  * ConcOut.N1_in_R1N1_Rab4_inact + parameters.kVegfR1N1Rab4a  * ConcOut.N1_in_R1N1_Rab4_ACT ...
                + parameters.kR2N1Rab4a  * ConcOut.N1_in_R2N1_Rab4_inact + parameters.kVegfR2N1Rab4a  * ConcOut.N1_in_R2N1_Rab4_ACT;
transferRab411  = parameters.kN1Rab4at11a  * (ConcOut.N1_Rab4  - ConcOut.N1_in_R1N1_Rab4_inact - ConcOut.N1_in_R1N1_Rab4_ACT - ConcOut.N1_in_R2N1_Rab4_inact - ConcOut.N1_in_R2N1_Rab4_ACT) ...   
                + parameters.kR1N1Rab4at11a  * ConcOut.N1_in_R1N1_Rab4_inact + parameters.kVegfR1N1Rab4at11a  * ConcOut.N1_in_R1N1_Rab4_ACT ...
                + parameters.kR2N1Rab4at11a  * ConcOut.N1_in_R2N1_Rab4_inact + parameters.kVegfR2N1Rab4at11a  * ConcOut.N1_in_R2N1_Rab4_ACT;
recyclingRab11  = parameters.kN1Rab11a  * (ConcOut.N1_Rab11  - ConcOut.N1_in_R1N1_Rab11_inact - ConcOut.N1_in_R1N1_Rab11_ACT - ConcOut.N1_in_R2N1_Rab11_inact - ConcOut.N1_in_R2N1_Rab11_ACT) ...   
                + parameters.kR1N1Rab11a  * ConcOut.N1_in_R1N1_Rab11_inact + parameters.kVegfR1N1Rab11a  * ConcOut.N1_in_R1N1_Rab11_ACT ...
                + parameters.kR2N1Rab11a  * ConcOut.N1_in_R2N1_Rab11_inact + parameters.kVegfR2N1Rab11a  * ConcOut.N1_in_R2N1_Rab11_ACT;
% NB would have to reformat the R1N1 terms if there are differential rate constants values for
% trafficking of VEGFR1 with PlGF bound vs VEGF bound (currently there are not)

N1rates.Surf_in_Synth(1:length(ConcOut.N1_surf)) = parameters.kN1prod;
% N1rates.Surf_in_Synth(2:length(ConcOut.N1_surf)) = parameters.kN1prod * perturbs.N1reductionByCHX;
N1rates.Surf_in_Rab4  = recyclingRab4;
N1rates.Surf_in_Rab11 = recyclingRab11;
N1rates.Surf_out_Rab4 = internalization;

N1rates.Rab4_in_Surf  = internalization;
N1rates.Rab4_out_Rab7 = degradation;
N1rates.Rab4_out_Surf = recyclingRab4;
N1rates.Rab4_out_Rab11= transferRab411;

N1rates.Rab11_in_Rab4 = transferRab411;
N1rates.Rab11_out_Surf= recyclingRab11;

N1rates.Surf_net = N1rates.Surf_in_Synth + N1rates.Surf_in_Rab4  + N1rates.Surf_in_Rab11 - N1rates.Surf_out_Rab4;
N1rates.Rab4_net = N1rates.Rab4_in_Surf  - N1rates.Rab4_out_Rab7 - N1rates.Rab4_out_Surf - N1rates.Rab4_out_Rab11;
N1rates.Rab11_net= N1rates.Rab11_in_Rab4 - N1rates.Rab11_out_Surf;

end
