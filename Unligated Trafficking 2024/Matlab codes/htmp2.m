function h1 = htmp2 (xsim,ysim,zsim,hmTitle,hmXlab,hmYlab,scalelims)

mapr = linspace(1, 1, 101);
mapg = linspace(1, 0, 101);
mapb = linspace(1, 0, 101);
map2  = [mapr' mapg' mapb'];

h1 = heatmap(xsim,ysim,zsim,'ColorLimits',scalelims); % (:,1,:) vegf and timepoints; % (1,:,:) plgf and timepoints
h1.Colormap = map2;
h1.FontSize = 14;
h1.Title = hmTitle;
h1.XLabel = hmXlab;
h1.YLabel = hmYlab;

end

