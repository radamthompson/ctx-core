function [ model, LPsoln ] = calcCoreFBA( model, expt )
%calcCoreFBA perform FBA on core C therm model
%   This function will take the core C. therm model in RAVEN format and set
%   a given experimental condition constraints.
%
%   model := Raven model
%   expt := experimental conditions
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

rxns = model.rxns;

iSUB = find(strcmp(rxns,'CEL1'));
iETOH = find(strcmp(rxns,'TRA1'));
iACE = find(strcmp(rxns,'TRA2'));
iLAC = find(strcmp(rxns,'TRA4'));
iFOR = find(strcmp(rxns,'TRA6'));
iBIO = find(strcmp(rxns,'BIO'));
iCO2 = find(strcmp(rxns,'TRA7'));
iH2 = find(strcmp(rxns,'TRA8'));
iNH3 = find(strcmp(rxns,'TRA9'));
iVAL = find(strcmp(rxns,'VAL5'));
iMAIN = find(strcmp(rxns,'MAINT'));

% Carbon balance
cCEL = 12; cGLC = 6; cETOH = 2; cACE = 2; cLAC = 3;
cFOR = 1; cSUCC = 4; cVAL = 5; cCO2 = 1; cALA = 3; cBIO = 37;

switch expt
    case 'wt'
        model=setParam(model,'ub',iSUB,3.74);
        model=setParam(model,'lb',iSUB,3.4158);
        
        model=setParam(model,'ub',iETOH,4.289);
        model=setParam(model,'lb',iETOH,4.081);
        %
        model=setParam(model,'ub',iACE,2.718);
        model=setParam(model,'lb',iACE,2.5465);
        %
        model=setParam(model,'ub',iLAC,0.193);
        model=setParam(model,'lb',iLAC,0.175);
        %
        model=setParam(model,'ub',iFOR,1.777);
        model=setParam(model,'lb',iFOR,1.760);
        %
        %       model=setParam(model,'ub',iH2,8.236);
        %       model=setParam(model,'lb',iH2,7.488);
        %
        %         model=setParam(model,'ub',iVAL,0.8963);
        %         model=setParam(model,'lb',iVAL,0.6560);
        
    case 'he'
        model=setParam(model,'ub','OPM7r',0);
        model=setParam(model,'lb','OPM7r',0);
        
        model=setParam(model,'ub','OPM8r',0);
        model=setParam(model,'lb','OPM8r',0);
        
        model=setParam(model,'ub','FEM14',0);
        model=setParam(model,'lb','FEM14',0);
        
        model=setParam(model,'ub',iSUB,3.94);
        model=setParam(model,'lb',iSUB,2.71);
        
        model=setParam(model,'ub',iETOH,9.316);
        model=setParam(model,'lb',iETOH,6.582);

        model=setParam(model,'ub',iACE,0.832);
        model=setParam(model,'lb',iACE,0.501);
        %
        model=setParam(model,'ub',iLAC,0.047);
        model=setParam(model,'lb',iLAC,0.042);
        %
        model=setParam(model,'ub',iFOR,3.108);
%        model=setParam(model,'lb',iFOR,2.921);
        %
        %       model=setParam(model,'ub',iH2,0.646);
        %       model=setParam(model,'lb',iH2,0.524);
        %
        %         model=setParam(model,'ub',iVAL,0.314);
        %         model=setParam(model,'lb',iVAL,0.309);
        
end



LPsoln = solveLP(model,3);