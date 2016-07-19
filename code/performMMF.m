function [ output ] = performMMF( model,obj,delrxn,fig_title )
% performMMF Perform single iteration of MMF algorithm
%   This function performs a rigidity analysis on a metatool model as a
%   single iteration of the Minimal Metabolic Functionality process (Trinh
%   et al., AEM, 2008).
%
%   model := Metatool metabolic model structure with species identifier
%   obj := Substrate preference; a string 'RXN_NAME'
%   delrxn := Cell array of reactions to delete on this iteration
%   fig_title := (Optional) If you want to print the figure to PDF, give it
%   a title
%
%   output := structure of rigidity stats
%
% R. Adam Thompson
% Trinh Lab
% University of Tennessee, Knoxville
% rthomp46@utk.edu
%
% Created: February 12, 2012 (First editions authored by Cong T. Trinh)
% Last edit: February 11, 2015


% Check for figure printing
if nargin > 3
    print_fig = 1;
else
    print_fig = 0;
end

% Make life a bit easier
ems = model.ems;
react_name = model.react_name;
species = model.species;
output = model;

% Only use modes with substrate of interest
if ~isempty(obj)
    for i = 1:length(obj)
        i_obj = find(strcmp(react_name,obj(i)));
        ems   = ems(:,find(ems(i_obj,:)~=0));
        if ~length(i_obj)
            fprintf('Error: Reaction %s is not in react_name\n',cell2mat(obj(i)))
            return
        end
    end
end

% Enforce deletion strategy
if ~isempty(delrxn)
    for i = 1:length(delrxn)
        i_delrxn = find(strcmp(react_name,delrxn(i)));
        ems   = ems(:,find(ems(i_delrxn,:)==0));
        if ~length(i_delrxn)
            fprintf('Error: Reaction %s is not in react_name\n',cell2mat(delrxn(i)))
            return
        end
    end
end

[m n] = size(ems);

if n == 0
    fprintf('There are no elementary modes with the given obj and delrxn\n')
    output.check = 0;
    return
else
    output.check = 1;
end

fprintf('Number of EMs: %6.0f\n',n)
output.ems = ems;

% Carbon balance
cCB = 12; cGLU = 6; cETOH = 2; cACE = 2; cLAC = 3;
cFOR = 1; cSUCC = 4; cVAL = 5; cCO2 = 1; cALA = 3;

% Get EM stats for this deletion set
switch species
    %%
    case 1
        % C. therm.
        iSUB = find(strcmp(react_name,'CEL1'));
        iETOH = find(strcmp(react_name,'TRA1'));
        iACE = find(strcmp(react_name,'TRA2'));
        iLAC = find(strcmp(react_name,'TRA4'));
        iFOR = find(strcmp(react_name,'TRA6'));
        iBIO = find(strcmp(react_name,'BIO'));
        iCO2 = find(strcmp(react_name,'TRA7'));
        iH2 = find(strcmp(react_name,'TRA8'));
        iVAL = find(strcmp(react_name,'VAL5'));
        
        output.ETOH_EMS = [];
        output.ACE_EMS = [];
        output.LAC_EMS = [];
        output.VAL_EMS = [];
        output.FOR_EMS = [];
        output.H2_EMS = [];
        output.BIO_EMS = [];
        output.EthBIO = [];
        
        cBIO = 37;
        
        % New yields
        output.YieldEtOH = ems(iETOH,:)*cETOH./(ems(iSUB,:)*cCB);
        output.YieldBIO = (ems(iBIO,:)*cBIO)./(ems(iSUB,:)*cCB);
        
        for i = 1:n
            if (ems(iETOH,i) ~= 0)
                output.ETOH_EMS(end+1) = i;
                if (ems(iBIO,i) ~= 0)
                    output.EthBIO(end+1) = i;
                end
            end
            if (ems(iACE,i) ~= 0)
                output.ACE_EMS(end+1) = i;
            end
            if (ems(iLAC,i) ~= 0)
                output.LAC_EMS(end+1) = i;
            end
            if (ems(iVAL,i) ~= 0)
                output.VAL_EMS(end+1) = i;
            end
            if (ems(iFOR,i) ~= 0)
                output.FOR_EMS(end+1) = i;
            end
            if (ems(iH2,i) ~= 0)
                output.H2_EMS(end+1) = i;
            end
            if (ems(iBIO,i) ~= 0)
                output.BIO_EMS(end+1) = i;
            end
            
        end
        
        % Initialize rigidity analysis
        output.range_ETOH = [];
        output.range_BIO = [];
        output.range_BOTH = [];
        output.modefrac = [];
        
        % Perform rigidity analysis
        for i = 1:length(react_name)
            ems_temp = ems;
            ems_temp = ems_temp(:, find(ems(i,:)==0));
            n_ems_temp = size(ems_temp,2);
            if ~isempty(find(ems_temp(iSUB,:)))&&~isempty(ems_temp)
                ndx = find(ems_temp(iSUB,:));
                YETOH = (ems_temp(iETOH,ndx)*cETOH)./(ems_temp(iSUB,ndx)*cCB);
                YBIO = (ems_temp(iBIO,ndx)*cBIO)./(ems_temp(iSUB,ndx)*cCB);
                YACE = (ems_temp(iACE,ndx)*cACE)./(ems_temp(iSUB,ndx)*cCB);
                output.range_ETOH(i,:) = [min(YETOH) max(YETOH)];
                output.range_BIO(i,:) = [min(YBIO) max(YBIO)];
                output.range_ACE(i,:) = [min(YACE) max(YACE)];
                [blah,eth] = find(YETOH);
                [blah,bio] = find(YBIO);
                BOTH = intersect(eth,bio);
                if ~isempty(BOTH)
                    output.range_BOTH(i,:) = [min(YETOH(BOTH)) max(YETOH(BOTH))];
                else
                    output.range_BOTH(i,:) = [0 0];
                end
            else
                output.range_ETOH(i,:) = [0 0];
                output.range_BIO(i,:) = [0 0];
                output.range_BOTH(i,:) = [0 0];
            end
            output.modefrac(i,1) = n_ems_temp/n;
        end
        
        % Sort new list by mode frac
        [mfrac ndx] = sort(output.modefrac,'ascend');
        ETOH = output.range_ETOH(ndx,:);
        BIO = output.range_BIO(ndx,:);
        rxn_name = react_name(ndx);
        iii = 1:1:length(react_name);
        max_sync = [];
        min_sync = [];
        bio_sync = [];
        
        for p = 1:length(mfrac)
            max_sync(p) = output.range_ETOH(ndx(p),2);
            min_sync(p) = output.range_ETOH(ndx(p),1);
            bio_sync(p) = output.range_BIO(ndx(p),2);
        end
        
        count = size(react_name);
        
        YieldBIO = output.YieldBIO;
        YieldEtOH = output.YieldEtOH;
        %%
        
    case 2
        % E. coli
        iSUB = find(strcmp(react_name,'GG1'));
        iETOH = find(strcmp(react_name,'TRA1'));
        iACE = find(strcmp(react_name,'TRA2'));
        iLAC = find(strcmp(react_name,'TRA4'));
        iFOR = find(strcmp(react_name,'TRA6'));
        iBIO = find(strcmp(react_name,'BIO'));
        iCO2 = find(strcmp(react_name,'TRA7'));
        iSUC = find(strcmp(react_name,'TRA5'));
        iO2 = find(strcmp(react_name,'OPM1'));
        iO2_2 = find(strcmp(react_name,'OPM2'));
        
        output.ETOH_EMS = [];
        output.ACE_EMS = [];
        output.LAC_EMS = [];
        output.SUC_EMS = [];
        output.FOR_EMS = [];
        output.ANAERO_EMS = [];
        output.BIO_EMS = [];
        output.EthBIO = [];
        
        cBIO = 45367;
        
        output.YieldEtOH = ems(iETOH,:)*cETOH./(ems(iSUB,:)*cGLU);
        output.YieldBIO = (ems(iBIO,:)*cBIO)./(ems(iSUB,:)*cGLU);
        
        for i = 1:n
            if (ems(iETOH,i) ~= 0)
                output.ETOH_EMS(end+1) = i;
                if (ems(iBIO,i) ~= 0)
                    output.EthBIO(end+1) = i;
                end
            end
            if (ems(iACE,i) ~= 0)
                output.ACE_EMS(end+1) = i;
            end
            if (ems(iLAC,i) ~= 0)
                output.LAC_EMS(end+1) = i;
            end
            if (ems(iSUC,i) ~= 0)
                output.SUC_EMS(end+1) = i;
            end
            if (ems(iFOR,i) ~= 0)
                output.FOR_EMS(end+1) = i;
            end
            if (ems(iBIO,i) ~= 0)
                output.BIO_EMS(end+1) = i;
            end
            if (ems(iO2,i) == 0 && ems(iO2_2,i) == 0)
                output.ANAERO_EMS(end+1) = i;
            end
            
        end
        
        % Initialize rigidity analysis
        output.range_ETOH = [];
        output.range_BIO = [];
        output.range_BOTH = [];
        output.modefrac = [];
        
        % Perform rigidity analysis
        for i = 1:length(react_name)
            ems_temp = ems;
            ems_temp = ems_temp(:, find(ems(i,:)==0));
            n_ems_temp = size(ems_temp,2);
            if ~isempty(find(ems_temp(iSUB,:)))&&~isempty(ems_temp)
                ndx = find(ems_temp(iSUB,:));
                YETOH = (ems_temp(iETOH,ndx)*cETOH)./(ems_temp(iSUB,ndx)*cGLU);
                YBIO = (ems_temp(iBIO,ndx)*cBIO)./(ems_temp(iSUB,ndx)*cGLU);
                output.range_ETOH(i,:) = [min(YETOH) max(YETOH)];
                output.range_BIO(i,:) = [min(YBIO) max(YBIO)];
                [blah,eth] = find(YETOH);
                [blah,bio] = find(YBIO);
                BOTH = intersect(eth,bio);
                if ~isempty(BOTH)
                    output.range_BOTH(i,:) = [min(BOTH) max(BOTH)];
                else
                    output.range_BOTH(i,:) = [0 0];
                end
                % fprintf('%i\n',i);
            else
                output.range_ETOH(i,:) = [0 0];
                output.range_BIO(i,:) = [0 0];
                output.range_BOTH(i,:) = [0 0];
            end
            output.modefrac(i,1) = n_ems_temp/n;
        end
        
        % Sort new list by mode frac
        [mfrac ndx] = sort(output.modefrac,'ascend');
        ETOH = output.range_ETOH(ndx,:);
        BIO = output.range_BIO(ndx,:);
        rxn_name = react_name(ndx);
        iii = 1:1:length(react_name);
        max_sync = [];
        min_sync = [];
        bio_sync = [];
        
        for p = 1:length(mfrac)
            max_sync(p) = output.range_ETOH(ndx(p),2);
            min_sync(p) = output.range_ETOH(ndx(p),1);
            bio_sync(p) = output.range_BIO(ndx(p),2);
        end
        
        count = size(react_name);
        
        YieldBIO = output.YieldBIO;
        YieldEtOH = output.YieldEtOH;
        %%
    case 3
        % Toy Network
        iSUB = find(strcmp(react_name,'R1'));
        iP = find(strcmp(react_name,'R4'));
        iD = find(strcmp(react_name,'R9'));
        
        output.YieldP = ems(iP,:)./ems(iSUB,:);
        output.YieldD = ems(iD,:)./ems(iSUB,:);
        
        output.P_EMS = find(output.YieldP);
        output.D_EMS = find(output.YieldD);
        
        % Initialize rigidity analysis
        output.range_P = [];
        output.range_D = [];
        output.modefrac = [];
        
        % Perform rigidity analysis
        for i = 1:length(react_name)
            ems_temp = ems;
            ems_temp = ems_temp(:, find(ems(i,:)==0));
            n_ems_temp = size(ems_temp,2);
            if ~isempty(find(ems_temp(iSUB,:)))&&~isempty(ems_temp)
                ndx = find(ems_temp(iSUB,:));
                Y_P = ems_temp(iP,ndx)./ems_temp(iSUB,ndx);
                Y_D = ems_temp(iD,ndx)./ems_temp(iSUB,ndx);
                output.range_P(i,:) = [min(Y_P) max(Y_P)];
                output.range_D(i,:) = [min(Y_D) max(Y_D)];
            else
                output.range_P(i,:) = [0 0];
                output.range_D(i,:) = [0 0];
            end
            output.modefrac(i,1) = n_ems_temp/n;
        end
        
        % Sort new list by mode frac
        [mfrac ndx] = sort(output.modefrac,'ascend');
        P = output.range_P(ndx,:);
        D = output.range_D(ndx,:);
        rxn_name = react_name(ndx);
        iii = 1:1:length(react_name);
        max_sync = [];
        min_sync = [];
        bio_sync = [];
        
        for p = 1:length(mfrac)
            max_sync(p) = output.range_P(ndx(p),2);
            min_sync(p) = output.range_P(ndx(p),1);
            bio_sync(p) = output.range_D(ndx(p),2);
        end
        
        count = size(react_name);
        
        YieldBIO = output.YieldD;
        YieldEtOH = output.YieldP;
        
        %%
end


%%
% Figure set up
fp = fillPage(gcf, 'margins', [0 0 0 0], 'papersize', [11 8.5]);
if print_fig == 1
    title(fig_title);
else
    title('MMF Results');
end

% Mode Fraction and Max Ethanol Yield per Reaction Deletion
subplot(2,1,1);
plot(iii, mfrac,'kx', iii, max_sync,'bo',iii,bio_sync,'g^',iii,min_sync,'r+');
xlabel('Reaction');switch species
    case 1
        xlim([1 60]); % C. therm.
        ylim([0 1]);
        ylabel('Y EtOH / Y BIO / Mode Fraction');
        
    case 2
        xlim([1 61]); % E. coli
        ylim([0 1]);
        ylabel('Y EtOH / Y BIO / Mode Fraction');
        
    case 3
        xlim([1 11]); % Toy Net
        ylim([0 2]);
        ylabel('Y_P / Y_D / Mode Fraction');
        
end
set(gca,'xtick',1:count);
set(gca,'xticklabel',rxn_name);
xticklabel_rotate([],90,[],'Fontsize',8);

% Biomass vs Ethanol Production
subplot(2,1,2);
plot(YieldBIO,YieldEtOH,'bx');
switch species
    case 1
        xlim([0 1]);
        ylim([0 1]);
        xlabel('Y_{BIO/CB} (C mol / C mol)');
        ylabel('Y_{EtOH/CB} (C mol / C mol)');
    case 2
        xlim([0 1]);
        ylim([0 1]);
        xlabel('Y_{BIO/GLU} (C mol / C mol)');
        ylabel('Y_{EtOH/GLU} (C mol / C mol)');
    case 3
        xlim([0 2]);
        ylim([0 2]);
        xlabel('Y_{D/A} (mol / mol)');
        ylabel('Y_{P/A} (mol / mol)');
end


if print_fig == 1
    print (gcf,'-dpdf','-r300',fig_title);
end


end

