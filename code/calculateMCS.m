function CNA_out = calculateMCS(model,minY,minRxn)
%calculateMCS Calculate MCS for network
% This function takes in a metabolic network and a minimum yield request
% then calculates the minimal cut set for the given inputs.  This version
% is designed to work with C. therm. networks with ethanol as the product
% of interest
%
%   model := Either a structure from output of performMMF or a string with
%   the name of the input metatool text file
%   minY := Specified starting minimum yield of desired product (ethanol),
%   optional, default 0.8
%   minRxn := Specified minimum number of reactions in the MCS, optional,
%   default 5
%
%   cutsets := output from 
%
% R. Adam Thompson
% Trinh Lab
% University of Tennessee, Knoxville
% rthomp46@utk.edu
%
% Original algorithm written 9/2012
% Last edit: March 9, 2015

% Check model input type and initialize. If string, calc EMs
if isa(model,'char')
    [cnap,errval] = CNAmetatool2MFNetwork(model);
    if errval == 0
        CNA_net = CNAsaveNetwork(cnap);
    else
        fprintf('****\nError in CNA model import\n*****\n')
    end
    [ems,rev,idx,ray] = CNAcomputeEFM(CNA_net,[],3,1,0,0,[],'All',{});
    CNA_net.ems = ems';
    CNA_net.react_name = strtrim(num2cell(CNA_net.reacID,2));
    CNA_net.species = 1;
    model = performMMF(CNA_net,{'CEL1'},[]);
end

% Simplify things
if nargin < 2
    minY = 0.8;
end
if nargin < 3
    minRxn = 5;
end

maxY = max(model.YieldEtOH);
yield = model.YieldEtOH;
allMax = yield == maxY;

react_name = model.react_name;
iCB = strcmp(react_name,'CEL1');
iETOH = strcmp(react_name,'TRA1');
iTAM = strcmp(react_name,'TAM1r');
iVAL1 = strcmp(react_name,'VAL1');
iVAL2 = strcmp(react_name,'VAL2r');
iVAL3 = strcmp(react_name,'VAL3r');
iVAL4 = strcmp(react_name,'VAL4r');
iVAL5 = strcmp(react_name,'VAL5');

% Set up MCS calculations
keep1 = yield >= (minY*maxY); 
keep2 = model.YieldBIO > 0;
keep3 = false(length(keep1),1);
for i = 1:length(keep1)
    if keep1(i) && keep2(i)
        keep3(i) = true;
    end
end
keep = allMax' + keep3;
keep = logical(keep);
toss = yield == 0;

goods = model.ems(:,keep)';
trash = model.ems(:,toss)';

sets2save.tabl2save = goods;
sets2save.min2save = 1;

%%%% Use below for when there are no ElMos needed to be kept
%CNA_out.cutsets = CNAcomputeCutsets(trash,'Inf',react_name,[],1);

%%%% Use below for when a set of ElMos should be kept
CNA_out.cutsets = CNAcomputeCutsets(trash,'Inf',react_name,sets2save,1);

% Prune MCSs
[m n] = size(CNA_out.cutsets);
cs = CNA_out.cutsets;
arr = [1:1:m]';
idx = [];
for i = 1:m
    if cs(i,iCB)
        idx = [idx; i];
    elseif cs(i,iETOH)
        idx = [idx; i];
    elseif cs(i,iTAM)
        idx = [idx; i];
    elseif cs(i,iVAL1)
        idx = [idx; i]; 
    elseif cs(i,iVAL2)
        idx = [idx; i];
    elseif cs(i,iVAL3)
        idx = [idx; i];
    elseif cs(i,iVAL4)
        idx = [idx; i];
    elseif cs(i,iVAL5)
        idx = [idx; i];
    end
end
idx = unique(idx);
prune = setdiff(arr,idx);

CNA_out.pruned = cs(prune,:);

% a = sum(CNA_test.pruned,1);
% a = num2str(a');
% [char(react_name), a]
