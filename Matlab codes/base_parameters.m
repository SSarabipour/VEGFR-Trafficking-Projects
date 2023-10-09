function parameters = base_parameters(model,baseparams)

if model == "Unligated_VEGFR_model_20230831"
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
        %Units; receptors/cell
        parameters.VEGFR1_0 = 18000; % total R1 calculated from 1800 rec/cell on surface/total fraction of 0.1
        parameters.VEGFR2_0 =  9607; % total R2 calculated from 4900 rec/cell on surface/total fraction of 0.51
        parameters.NRP1_0   = 91891; % total R1 calculated from 68000 rec/cell on surface/total fraction of 0.74

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
        kN1R1on = 8e-4; %1.66e-2; % units = (receptors/um2)-1.s-1
        parameters.kN1R1on_surf  = kN1R1on/SA_surface;
        parameters.kN1R1on_rab4  = kN1R1on/SA_rab4;
        parameters.kN1R1on_rab11 = kN1R1on/SA_rab11;
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
