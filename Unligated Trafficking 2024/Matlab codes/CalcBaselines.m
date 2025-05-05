function baselines = CalcBaselines(ConcOut)

outnames = fieldnames(ConcOut);

for j = 1:length(outnames)
	baselines.(outnames{j}) = ConcOut.(outnames{j})(end);
end

end
