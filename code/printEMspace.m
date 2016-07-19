function [ output ] = printEMspace( ems, react_name, fig_title )
%plotEMspace 
%   This function will take a set of EMs, calculate the yield of product
%   per substrate, then plot Yield_{P/S} vs Yield_{Bio/S} and export to pdf
%
%   Set up for C therm. core network (cellobiose)
%
%   Inputs:
%       ems - A set of ems for the network in question
%       react_name - A cell array with the reaction names in the network
%       fig_title - A string which will be the name of the pdf produced
%
%
% R. Adam Thompson
% Trinh Lab
% University of Tennessee, Knoxville
% rthomp46@utk.edu

cCB = 12;    
cETOH = 2;
cBIO = 4;

iCEL = find(strcmp(react_name,'CEL1'));
iBIO = find(strcmp(react_name,'BIO'));
iETOH = find(strcmp(react_name,'TRA1'));

output.YieldEtOH = ems(iETOH,:)*cETOH./(ems(iCEL,:)*cCB);
output.YieldBIO = (ems(iBIO,:)*cBIO)./(ems(iCEL,:)*cCB);

output.YieldEtOH(isnan(output.YieldEtOH))=0;
output.YieldBIO(isnan(output.YieldBIO))=0;


figure
plot(output.YieldBIO,output.YieldEtOH,'b*');
xlabel('Y_{BIO/CB} (C mol / C mol)');
ylabel('Y_{EtOH/CB} (C mol / C mol)');
xlim([0 1]);
ylim([0 1]);

set(gca,'XTick',[0:0.1:1.0],'YTick',[0:0.1:1.0])
set(gca,'FontSize',14)
set(findall(gcf,'type','text'),'fontSize',16)
fp = fillPage(gcf, 'margins', [0 0 0 0], 'papersize', [8.5 8.5]);
print (gcf,'-dpdf','-r300',fig_title);

end

