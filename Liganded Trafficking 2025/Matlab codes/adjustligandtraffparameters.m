function parameters = adjustligandtraffparameters(parameters,scanparams)


%% 4. TRAFFICKING PARAMETERS 
% (all first order processes, all s-1)

%% INTERNALIZATION
 
parameters.kVegfR1Rab5a = parameters.kR1Rab5a * scanparams.adjkintR1vegf;
parameters.kVegfR1N1Rab5a = parameters.kVegfR1Rab5a;
parameters.kPlgfR1Rab5a = parameters.kR1Rab5a * scanparams.adjkintR1plgf;
parameters.kPlgfR1N1Rab5a = parameters.kPlgfR1Rab5a;
parameters.kVegfR2Rab5a = parameters.kR2Rab5a * 3.0 * scanparams.adjkintR2vegf;
parameters.kVegfR2N1Rab5a = parameters.kVegfR2Rab5a;
 
parameters.kVegfN1Rab5a = parameters.kN1Rab5a * scanparams.adjkintN1vegf;
parameters.kPlgfN1Rab5a = parameters.kN1Rab5a * scanparams.adjkintN1plgf;
 

%% RECYCLING via Rab4
 
parameters.kVegfR1Rab4a = parameters.kR1Rab4a * scanparams.adjkrec4R1vegf;
parameters.kVegfR1N1Rab4a = parameters.kVegfR1Rab4a;
parameters.kPlgfR1Rab4a = parameters.kR1Rab4a * scanparams.adjkrec4R1plgf;
parameters.kPlgfR1N1Rab4a = parameters.kPlgfR1Rab4a;
parameters.kVegfR2Rab4a = parameters.kR2Rab4a * scanparams.adjkrec4R2vegf;
parameters.kVegfR2N1Rab4a = parameters.kVegfR2Rab4a;
 
parameters.kVegfN1Rab4a = parameters.kN1Rab4a * scanparams.adjkrec4N1vegf;
parameters.kPlgfN1Rab4a = parameters.kN1Rab4a * scanparams.adjkrec4N1plgf;
 
%% TRANSFER Rab4=>Rab11
 
parameters.kVegfR1Rab4at11a = parameters.kR1Rab4at11a * scanparams.adjkrec4t11R1vegf;
parameters.kVegfR1N1Rab4at11a = parameters.kVegfR1Rab4at11a;
parameters.kPlgfR1Rab4at11a = parameters.kR1Rab4at11a * scanparams.adjkrec4t11R1plgf;
parameters.kPlgfR1N1Rab4at11a = parameters.kPlgfR1Rab4at11a;
parameters.kVegfR2Rab4at11a = parameters.kR2Rab4at11a * scanparams.adjkrec4t11R2vegf;
parameters.kVegfR2N1Rab4at11a = parameters.kVegfR2Rab4at11a;
 
parameters.kVegfN1Rab4at11a = parameters.kN1Rab4at11a * scanparams.adjkrec4t11N1vegf;
parameters.kPlgfN1Rab4at11a = parameters.kN1Rab4at11a * scanparams.adjkrec4t11N1plgf;
 
%% RECYCLING via Rab11
 
parameters.kVegfR1Rab11a = parameters.kR1Rab11a * scanparams.adjkrec11R1vegf;
parameters.kVegfR1N1Rab11a = parameters.kVegfR1Rab11a;
parameters.kPlgfR1Rab11a = parameters.kR1Rab11a * scanparams.adjkrec11R1plgf;
parameters.kPlgfR1N1Rab11a = parameters.kPlgfR1Rab11a;
parameters.kVegfR2Rab11a = parameters.kR2Rab11a * scanparams.adjkrec11R2vegf;
parameters.kVegfR2N1Rab11a = parameters.kVegfR2Rab11a;
 
parameters.kVegfN1Rab11a = parameters.kN1Rab11a * scanparams.adjkrec11N1vegf;
parameters.kPlgfN1Rab11a = parameters.kN1Rab11a * scanparams.adjkrec11N1plgf;
 
%% DEGRADATION
 
parameters.kVegfR1Rab4at7a = parameters.kR1Rab4at7a * scanparams.adjkdegR1vegf;
parameters.kVegfR1N1Rab4at7a = parameters.kVegfR1Rab4at7a;
parameters.kPlgfR1Rab4at7a = parameters.kR1Rab4at7a * scanparams.adjkdegR1plgf;
parameters.kPlgfR1N1Rab4at7a = parameters.kPlgfR1Rab4at7a;
parameters.kVegfR2Rab4at7a = parameters.kR2Rab4at7a * scanparams.adjkdegR2vegf;
parameters.kVegfR2N1Rab4at7a = parameters.kVegfR2Rab4at7a;
 
parameters.kVegfN1Rab4at7a = parameters.kN1Rab4at7a * scanparams.adjkdegN1vegf;
parameters.kPlgfN1Rab4at7a = parameters.kN1Rab4at7a * scanparams.adjkdegN1plgf;
 
end
