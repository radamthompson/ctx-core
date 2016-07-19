function [ output ] = METAFoR( em_struc, flux_struc, fig_title )
%METAFoR Perform Metabolic Flux Ratio Analysis
%   This function will take a flux structure (output from getFluxDistr) and
%   perform Metabolic Flux Ratio Analysis on the flux distribution vector.
%   For details on the method, see Sauer et al. J. Bacteriol., 1999,
%   6679-6688.
%
%   This version is for the C. thermocellum network, and so the reactions
%   of interest differ from the Sauer reference.
%
%       em_space := Output structure from runMetatool or performMMF
%       flux_struc := Output structure from getFluxDistr
%       fig_title := String containing the name of the output pdf
%       (optional)
%
%       output := METAFoR statistics structure
%
%
% R. Adam Thompson
% Trinh Lab
% University of Tennessee, Knoxville
% rthomp46@utk.edu
%
% Created: March 8, 2015
% Last edit: 


% Check for figure printing
if nargin > 1
    print_fig = 1;
else
    print_fig = 0;
end

% Initialize
react_name = em_struc.react_name;
r_orig = flux_struc.r;

tPYR = 4; % Theoretical max from 1 CB (mol / mol)


% Get indices for reactions of interest
iCB = find(strcmp(react_name,'CEL1'));
iOAADC = find(strcmp(react_name,'ANA4'));
iMAE = find(strcmp(react_name,'ANA3'));
iPPDK = find(strcmp(react_name,'GG11r'));
iECH = find(strcmp(react_name,'FEM14'));
iBIF = find(strcmp(react_name,'OPM8r'));
iFEH2 = find(strcmp(react_name,'OPM7r'));
iH2 = find(strcmp(react_name,'TRA8'));
iPFOR = find(strcmp(react_name,'GG13'));
iPFL = find(strcmp(react_name,'FEM1'));
iACE = find(strcmp(react_name,'FEM7'));
iETH = find(strcmp(react_name,'FEM5'));
iTCA = find(strcmp(react_name,'TCA1'));
iLAC = find(strcmp(react_name,'FEM3'));
iRNF = find(strcmp(react_name,'OPM4r'));
iNFN = find(strcmp(react_name,'OPM9r'));
    iETH2 = find(strcmp(react_name,'FEM15')); % NADPH dependent reaction
if iETH2
    eth2 = 1;
else
    eth2 = 0;
end

% Normalize r to cellobiose uptake
r = r_orig./r_orig(iCB);
%r = abs(r);

% Get total fluxes
rPYR = r(iOAADC) + r(iMAE) + r(iPPDK);
rH2 = r(iECH) + r(iBIF) - r(iFEH2);
rFDRD = r(iECH) + r(iRNF) + r(iBIF) + r(iNFN);
rACoA = r(iPFL) + r(iPFOR);
if eth2 == 1
rFERM = r(iACE) + r(iETH) + r(iETH2) + r(iTCA);
else
    rFERM = r(iACE) + r(iETH) + r(iTCA);
end

% Calculate flux ratios
fMALS = (r(iOAADC) + r(iMAE))/rPYR;
fPPDK = r(iPPDK) / rPYR;

fECH = r(iECH)/rFDRD;
fRNF = r(iRNF)/rFDRD;
fBIF = r(iBIF)/rFDRD;
fNFN = r(iNFN)/rFDRD;

fECH_2 = r(iECH)/rH2;
fBIF_2 = r(iBIF)/rH2;
fFEH2 = r(iFEH2)/rH2;

fPFL = r(iPFL)/rACoA;
fPFOR = r(iPFOR)/rACoA;

fLAC = r(iLAC)/rPYR;

fTCA = r(iTCA)/rFERM;

fACE = r(iACE)/rFERM;
if eth2 == 1
fETH = (r(iETH)+r(iETH2))/rFERM;
else
    fETH = r(iETH)/rFERM;
end

% Make output and printer friendly
name_vec = {'Mal Shunt';'PPDK';'ECH';'RNF';'BIF';'NFN';'ECH';'BIF';'Fe-H2';'PFL';'PFOR';'LDH';'TCA Cycle';'PTA-ACK';'AdhE'};
f_vec = [fMALS; fPPDK; fECH; fRNF; fBIF; fNFN; fECH_2; fBIF_2; fFEH2; fPFL; fPFOR; fLAC; fTCA; fACE; fETH];

output.name_vec = name_vec;
output.f_vec = f_vec;

output.r_orig = flux_struc.r;

end

