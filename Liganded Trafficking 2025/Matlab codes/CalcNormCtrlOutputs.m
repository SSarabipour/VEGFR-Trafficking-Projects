function ConcOut=CalcNormCtrlOutputs(ConcIn,ConcCtrl)

ConcOut = ConcIn;
outnames = fieldnames(ConcIn);
for i = 1:length(outnames)
    ConcOut.(outnames{i}) = ((ConcIn.(outnames{i}))./ConcCtrl.(outnames{i}))*100; 
    % divide each timepoint by the equivalent timepoint in the contrl (
    % (and *100 to get percentage of control)
end
end





