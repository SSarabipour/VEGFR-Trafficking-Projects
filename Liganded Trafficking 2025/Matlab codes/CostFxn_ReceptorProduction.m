function CostOut = CostFxn_ReceptorProduction(p2,targets,parameters,model)
%
% INPUTS: p2 = current guesses for receptor production rates
%         parameters = parameter values (other than production rates)
%         targets = targets for cell-surface VEGFR1, VEGFR2, NRP1 totals
%
% OUTPUTS: CostOut = three-element array of errors (i.e. differences between 
%               the simulated and target VEGFR1, VEGFR2, NRP1 surface expressions)
%

p2(targets==0) = 0; 

parameters.kR1prod = p2(1);
parameters.kR2prod = p2(2);
parameters.kN1prod = p2(3);

[timepoints, species_out, observables_out] = SimToSteadyState(parameters,model);

ConcOut = CalcOutputsL(observables_out);

observSim(1) = ConcOut.R1_surf(end); 
observSim(2) = ConcOut.R2_surf(end); 
observSim(3) = ConcOut.N1_surf(end); 

for j=1:3
	CostOut(j) = ((observSim(j) - targets(j)));
end

% CostOut

end
