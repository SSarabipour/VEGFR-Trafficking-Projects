function parameters = base_parameters(model,baseparams)

if model == "Unligated_VEGFR_model_20230321"
    runmodel = 1
end

switch runmodel
    case 1

        % set up the parameter structure (using a consistent order)
        a = fileread("Paramlist.txt");
        pn = strsplit(a);
        for i=1:length(pn)
            parameters.(pn{i}) = 0;
        end

        %% IMPORTANT: BASE UNITS OF SIMULATION - molecules per um2 of cell surface (typically 1000 um2/cell)

        %% INITIAL CONDITIONS

        %% Initial receptor levels
        %Units; receptors/um^2
        parameters.VEGFR1_0 = 18; % total R1 calculated from 1800 rec/cell on surface/total fraction of 0.1
        parameters.VEGFR2_0 =  9.607; % total R2 calculated from 4900 rec/cell on surface/total fraction of 0.51
        parameters.NRP1_0   = 91.891; % total N1 calculated from 68000 rec/cell on surface/total fraction of 0.74

        %% Receptor-Receptor association/dissociation parameters 
        % Units: kon = (receptors/um2)-1.s-1; koff = s-1
        parameters.kR1R1on  = 0.008;
        parameters.kR1R1off = 0.01;
        parameters.kR2R2on  = 0.002;
        parameters.kR2R2off = 0.01;
        parameters.kN1R1on  = 0.0166;
        parameters.kN1R1off = 0.01;

        %% TRAFFICKING PARAMETERS (all first order processes, all s-1)

        p1 = readmatrix(baseparams); % import unliganded parameters from csv file

        %% INTERNALIZATION
            parameters.kR1Rab5a = p1(1);
            parameters.kR2Rab5a = p1(6);
            parameters.kN1Rab5a = p1(11);
            parameters.kR1N1Rab5a = p1(1); % we set this rate the same as R1 internalization

        %% DEGRADATION
            parameters.kR1Rab4at7a = p1(2);
            parameters.kR2Rab4at7a = p1(7);
            parameters.kN1Rab4at7a = p1(12);
            parameters.kR1N1Rab4at7a = p1(2); % we set this rate the same as R1 degradation

        %% RECYCLING via Rab4
            parameters.kR1Rab4a = p1(3);
            parameters.kR2Rab4a = p1(8);
            parameters.kN1Rab4a = p1(13);
            parameters.kR1N1Rab4a = p1(3); % we set this rate the same as R1 recycling

        %% TRANSFER Rab4=>Rab11
            parameters.kR1Rab4at11a = p1(4);
            parameters.kR2Rab4at11a = p1(9);
            parameters.kN1Rab4at11a = p1(14);
            parameters.kR1N1Rab4at11a = p1(4); % we set this rate the same as R1 transfer

        %% RECYCLING via Rab11
            parameters.kR1Rab11a = p1(5);
            parameters.kR2Rab11a = p1(10);
            parameters.kN1Rab11a = p1(15);
            parameters.kR1N1Rab11a = p1(5); % we set this rate the same as R1 recycling

        %% RECEPTOR PRODUCTION
            parameters.kR1prod = p1(16);
            parameters.kR2prod = p1(17);
            parameters.kN1prod = p1(18);
 
end
end
