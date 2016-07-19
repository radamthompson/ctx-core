function [ ETH,BIO ] = findHullAndPlot( MetatoolModel, fig_title )
%findHullAndPlot
%   This function will take a data structure from the output of
%   performMMF and output the convex hull of the phenotypic space
%   between ethanol and biomass yield. In addition, it can optionally
%   output a figure of the entire Eth/Bio yield space
%
% R. Adam Thompson
% Trinh Lab
% University of Tennessee, Knoxville
% rthomp46@utk.edu

printFig = 1;
if nargin==1
    printFig = 0;
end

ETH = MetatoolModel.YieldEtOH;
BIO = MetatoolModel.YieldBIO;
hull = convhull(BIO,ETH);
ETH = MetatoolModel.YieldEtOH(hull)';
BIO = MetatoolModel.YieldBIO(hull)';

figure
plot(MetatoolModel.YieldBIO,MetatoolModel.YieldEtOH,'b*');
xlabel('Y_{BIO/CB} (C mol / C mol)');
ylabel('Y_{EtOH/CB} (C mol / C mol)');
xlim([0 1]);
ylim([0 1]);

set(gca,'XTick',[0:0.1:1.0],'YTick',[0:0.1:1.0])
set(gca,'FontSize',14)
set(findall(gcf,'type','text'),'fontSize',16)
fp = fillPage(gcf, 'margins', [0 0 0 0], 'papersize', [8.5 8.5]);

if printFig == 1
    print (gcf,'-dpdf','-r300',fig_title);
end

end

