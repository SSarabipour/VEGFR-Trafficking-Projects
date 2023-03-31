function h1 = htmp (xsim,ysim,zsim,hmTitle,hmXlab,hmYlab)

mapr = [linspace(0.33,1,101)  linspace(.99,0,100)];
mapg = [linspace(0,1,101)   linspace(.99,.33,100)];
mapb = [linspace(0.33,1,101) linspace(0.99,0,100)];
map  = [mapr' mapg' mapb'];

h1 = heatmap(xsim,ysim,zsim); % (:,1,:) vegf and timepoints; % (1,:,:) plgf and timepoints
h1.Colormap = map;
h1.FontSize = 14;
h1.Title = hmTitle;
h1.XLabel = hmXlab;
h1.YLabel = hmYlab;

end

