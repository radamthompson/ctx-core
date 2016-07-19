function [ output ] = plotEMspace( ems, react_name )
%plotEMspace 
%   This function will take a set of EMs, calculate the yield of product
%   per substrate, then plot Yield_{P/S} vs Yield_{Bio/S}

cCB = 12;    
cETOH = 2;
cBIO = 4;

iCEL = find(strcmp(react_name,'CEL1'));
iBIO = find(strcmp(react_name,'BIO'));
iETOH = find(strcmp(react_name,'TRA1'));

output.YieldEtOH = ems(iETOH,:)*cETOH./(ems(iCEL,:)*cCB);
output.YieldBIO = (ems(iBIO,:)*cBIO)./(ems(iCEL,:)*cCB);

figure
plot(output.YieldBIO,output.YieldEtOH,'b*');
xlabel('Y_{BIO/CB} (C mol / C mol)');
ylabel('Y_{EtOH/CB} (C mol / C mol)');
xlim([0 1]);
ylim([0 1]);

end

