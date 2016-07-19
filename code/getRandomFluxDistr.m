function [ model,FluxSS ] = getRandomFluxDistr( model, r_min, r_max, species )
%getRandomFluxDistr Perform Metabolic Flux Analysis
%   This function will take a Metatool model structure and a set of
%   measured flux bounds to find a random flux distribution vector within
%   experimental bounds as part of the creation of experimentally bound
%   ensemble. The ensemble is used in the analysis of representative flux
%   vectors. The algorithm currently uses the Moore-Penrose pseudo-inverse
%   method as described in Poolman et al.
%
%   model := Metatool model structure, requires elementary modes to be in
%   field model.ems
%   r_min := Array of min measured flux values, must be in the same order as
%   model.ext_met and consistent in supplied values as r_max.
%   r_max := Array of max measured flux values, must be in the same order as
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
% Created: October 6, 2015 
% Last edit: 

react_name = model.react_name;
model.species = species;

% Remove empty flux values
r_bool = ~cellfun(@isempty,r_min);
r_min = r_min(r_bool);
r_max = r_max(r_bool);

% Generate random r_m
n = length(r_min);
rand_arr = rand(n,1);
r_m = cell(n,1);

for i = 1:n
    vari = r_max{i} - r_min{i};
    r_m{i} = r_max{i} - (vari * rand_arr(i));
end

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