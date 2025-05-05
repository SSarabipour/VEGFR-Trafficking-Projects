function ConcOut=CalcNormOutputs(ConcIn);

ConcOut = ConcIn;
outnames = fieldnames(ConcIn);
for i = 1:length(outnames)
    ConcOut.(outnames{i}) = ((ConcIn.(outnames{i}))./ConcIn.(outnames{i})(1))*100; 
    % divide each timepoint by initial timepoint as the control and *100 to
    % get percentage of control
end
end





