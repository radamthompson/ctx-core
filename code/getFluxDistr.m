function [ model,FluxSS ] = getFluxDistr( model, r_m, species )
%getFluxDistr Perform Metabolic Flux Analysis
%   This function will take a Metatool model structure and a set of
%   measured fluxes to find a flux distribution vector using MFA. The
%   algorithm currently uses the Moore-Penrose pseudo-inverse method.
%
%   model := Metatool model structure, requires elementary modes to be in
%   field model.ems
%   r_m := Array of measured flux values, must be in the same order as
%   model.ext_met
%   species := integer to discern the network being analyzed
%           1 := C. thermocellum
%           2 := E. coli
%           3 := Toy network
%
%   FluxSS := Structure with flux distribution calculated, yields, and EMs for specific
%   metabolites
%
% R. Adam Thompson
% Trinh Lab
% University of Tennessee, Knoxville
% rthomp46@utk.edu
%
% Created: September 28, 2012 (First editions authored by Cong T. Trinh)
% Last edit: February 11, 2015

react_name = model.react_name;
model.species = species;

% Remove empty flux values
r_bool = ~cellfun(@isempty,r_m);
r_m_in = r_m;
r_m = r_m(r_bool);

% Carbon balance
cCEL = 12; cGLC = 6; cETOH = 2; cACE = 2; cLAC = 3;
cFOR = 1; cSUCC = 4; cVAL = 5; cCO2 = 1; cALA = 3;


% Get indicies of transport reactions
i_rm = zeros(length(model.ext_met),1);

for i=1:length(model.ext_met)
    idx = find(model.ext(i,:));
    if length(idx) > 1
        fprintf('Caution, multiple reactions use metabolite %s\n',model.ext_met{i})
        i_rm(i) = idx(1);
    else
        i_rm(i) = idx;
    end
end

% Match rxns to measured fluxes
i_rm_full = i_rm;
i_rm = i_rm(r_bool);
ext_met = model.ext_met(r_bool);

% Set up indicies for specific species/networks
switch species
    case 1
        
        iSUB = i_rm_full(find(strcmp(ext_met,'CB_ext')));
        iETOH = i_rm_full(find(strcmp(ext_met,'ETOH_ext')));
        iACE = i_rm_full(find(strcmp(ext_met,'ACE_ext')));
        iLAC = i_rm_full(find(strcmp(ext_met,'LAC_ext')));
        iFOR = i_rm_full(find(strcmp(ext_met,'FOR_ext')));
        iBIO = i_rm_full(find(strcmp(ext_met,'BIOMASS')));
        iCO2 = i_rm_full(find(strcmp(ext_met,'CO2_ext')));
        iH2 = i_rm_full(find(strcmp(ext_met,'H2_ext')));
        iNH3 = i_rm_full(find(strcmp(ext_met,'NH3_ext')));
        iVAL = i_rm_full(find(strcmp(ext_met,'VAL_ext')));
        iMAIN = i_rm_full(find(strcmp(ext_met,'ATP_maint')));
        
        FluxSS.i_rm = i_rm;
        cBIO = 37; % For C. therm. network
        
        ndx = find(model.ems(iSUB,:)~=0);
        ems_new = model.ems(:,ndx);
    case 2
        iSUB = i_rm_full(find(strcmp(model.ext_met,'GLC_ext')));
        iETOH = i_rm_full(find(strcmp(model.ext_met,'ETOH_ext')));
        iACE = i_rm_full(find(strcmp(model.ext_met,'ACE_ext')));
        iLAC = i_rm_full(find(strcmp(model.ext_met,'LAC_ext')));
        iFOR = i_rm_full(find(strcmp(model.ext_met,'FOR_ext')));
        iBIO = i_rm_full(find(strcmp(model.ext_met,'BIOMASS')));
        iCO2 = i_rm_full(find(strcmp(model.ext_met,'CO2_ext')));
        iH2 = i_rm_full(find(strcmp(model.ext_met,'H2_ext')));
        iNH3 = i_rm_full(find(strcmp(model.ext_met,'NH3_ext')));
        iSUCC = i_rm_full(find(strcmp(model.ext_met,'SUCC_ext')));
        iMAIN = i_rm_full(find(strcmp(model.ext_met,'ATP_main')));
        
        FluxSS.i_rm = i_rm;
        cBIO = 45367; % For E. coli
        
        ndx = find(model.ems(iSUB,:)~=0);
        ems_new = model.ems(:,ndx);
        
    case 3
        % Toy Network
        iSUB = find(strcmp(react_name,'R1'));
        iP = find(strcmp(react_name,'R4'));
        iB = find(strcmp(react_name,'R8r'));
        iD = find(strcmp(react_name,'R9'));
        
        FluxSS.i_rm = i_rm;
        
        ndx = find(model.ems(iSUB,:)~=0);
        ems_new = model.ems(:,ndx);
end

switch species
    case 1
        % Overall Yield Calculations, uses all EMs, even if product is 0
        FluxSS.Yield.ETOH = model.ems(iETOH,:)*cETOH./(model.ems(iSUB,:)*cCEL);
        FluxSS.max.ETOH = max(FluxSS.Yield.ETOH);
        
        FluxSS.Yield.ACE = model.ems(iACE,:)*cACE./(model.ems(iSUB,:)*cCEL);
        FluxSS.max.ACE = max(FluxSS.Yield.ACE);
        
        FluxSS.Yield.LAC = model.ems(iLAC,:)*cLAC./(model.ems(iSUB,:)*cCEL);
        FluxSS.max.LAC = max(FluxSS.Yield.LAC);
        
        FluxSS.Yield.FOR = model.ems(iFOR,:)*cFOR./(model.ems(iSUB,:)*cCEL);
        FluxSS.max.FOR = max(FluxSS.Yield.FOR);
        
        FluxSS.Yield.H2 = model.ems(iH2,:)./(model.ems(iSUB,:));
        FluxSS.max.H2 = max(FluxSS.Yield.H2);
        
        FluxSS.Yield.BIO = model.ems(iBIO,:)*cBIO./(model.ems(iSUB,:)*cCEL);
        FluxSS.max.BIO = max(FluxSS.Yield.BIO);
        
        %if exist(iVAL)
            FluxSS.Yield.VAL = model.ems(iVAL,:)*cVAL./(model.ems(iSUB,:)*cCEL);
            FluxSS.max.VAL = max(FluxSS.Yield.VAL);
       %end
        
        FluxSS.Yield.CO2 = model.ems(iCO2,:)*cCO2./(model.ems(iSUB,:)*cCEL);
        FluxSS.max.CO2 = max(FluxSS.Yield.CO2);
        
        FluxSS.ETOH_EMS = find(FluxSS.Yield.ETOH);
        FluxSS.ACE_EMS = find(FluxSS.Yield.ACE);
        FluxSS.LAC_EMS = find(FluxSS.Yield.LAC);
        FluxSS.FOR_EMS = find(FluxSS.Yield.FOR);
        FluxSS.H2_EMS = find(FluxSS.Yield.H2);
        FluxSS.BIO_EMS = find(FluxSS.Yield.BIO);
        %if exist(iVAL)
            FluxSS.VAL_EMS = find(FluxSS.Yield.VAL);
        %end
        FluxSS.CO2_EMS = find(FluxSS.Yield.CO2);
        
    case 2
        
        FluxSS.Yield.ETOH = model.ems(iETOH,:)*cETOH./(model.ems(iSUB,:)*cGLC);
        FluxSS.max.ETOH = max(FluxSS.Yield.ETOH);
        
        FluxSS.Yield.ACE = model.ems(iACE,:)*cACE./(model.ems(iSUB,:)*cGLC);
        FluxSS.max.ACE = max(FluxSS.Yield.ACE);
        
        FluxSS.Yield.LAC = model.ems(iLAC,:)*cLAC./(model.ems(iSUB,:)*cGLC);
        FluxSS.max.LAC = max(FluxSS.Yield.LAC);
        
        FluxSS.Yield.FOR = model.ems(iFOR,:)*cFOR./(model.ems(iSUB,:)*cGLC);
        FluxSS.max.FOR = max(FluxSS.Yield.FOR);
        
        FluxSS.Yield.H2 = model.ems(iH2,:)./(model.ems(iSUB,:));
        FluxSS.max.H2 = max(FluxSS.Yield.H2);
        
        FluxSS.Yield.BIO = model.ems(iBIO,:)*cBIO./(model.ems(iSUB,:)*cGLC);
        FluxSS.max.BIO = max(FluxSS.Yield.BIO);
        
        FluxSS.Yield.CO2 = model.ems(iCO2,:)*cCO2./(model.ems(iSUB,:)*cGLC);
        FluxSS.max.CO2 = max(FluxSS.Yield.CO2);
        
        FluxSS.Yield.SUCC = model.ems(iSUCC,:)*cSUCC./(model.ems(iSUB,:)*cGLC);
        FluxSS.max.SUCC = max(FluxSS.Yield.SUCC);
        
        FluxSS.ETOH_EMS = find(FluxSS.Yield.ETOH);
        FluxSS.ACE_EMS = find(FluxSS.Yield.ACE);
        FluxSS.LAC_EMS = find(FluxSS.Yield.LAC);
        FluxSS.SUCC_EMS = find(FluxSS.Yield.SUCC);
        FluxSS.FOR_EMS = find(FluxSS.Yield.FOR);
        FluxSS.H2_EMS = find(FluxSS.Yield.H2);
        FluxSS.BIO_EMS = find(FluxSS.Yield.BIO);
        FluxSS.CO2_EMS = find(FluxSS.Yield.CO2);
        
    case 3
        FluxSS.Yield.P = model.ems(iP,:)./model.ems(iSUB,:);
        FluxSS.Yield.D = model.ems(iD,:)./model.ems(iSUB,:);
        
        FluxSS.max.P = max(FluxSS.Yield.P);
        FluxSS.max.D = max(FluxSS.Yield.D);
        
        FluxSS.P_ems = find(FluxSS.Yield.P);
        FluxSS.D_ems = find(FluxSS.Yield.D);
        
end

% Normalize EMs to substrate uptake rate
for i = 1:size(ems_new,2)
    r = cell2mat(r_m(iSUB));
    ems_new(:,i) = ems_new(:,i)*r/ems_new(iSUB,i);
end

EM_m = ems_new(FluxSS.i_rm,:);
FluxSS.rm = cell2mat(r_m);

%Perform flux distribution calculations - CTT
for i = 1:size(EM_m,2)
    w = pinv(EM_m)*FluxSS.rm;
    negative = find(w<0);
    if isempty(negative)
        break
    else
        EM_m(:,negative) = [];
        ems_new(:,negative) = [];
        ndx(negative) = [];
    end
end

FluxSS.w = w;
FluxSS.ndx_ems = ndx;
FluxSS.r = ems_new*FluxSS.w;
FluxSS.err = norm(FluxSS.r(FluxSS.i_rm)-FluxSS.rm,2)/norm(FluxSS.rm)*100;
FluxSS.calc_rm = FluxSS.r(FluxSS.i_rm);



end

