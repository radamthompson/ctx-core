%% yieldParse
% Adam Thompson
% 5/14/2014
% Trinh Lab, University of Tennessee

function [output,iEMs,newEMs] = yieldParse(input, EMs)
%%
% This script is designed to take in a vector of yields for a product of interest given the metabolic
% network. This can be for any product. The format works to use the output
% of the Ctherm_flux_calc.m script.
%
% output := an n x 1 vector corresponding to v_rep of high yielding EMs
% normalized to cellobiose uptake. 
% iEMs := indices of high-yielding elementary modes
%
% input := 1 x q vector of yields for the product of interest where q is the number of elementary modes
% in the network
%
% EMs := m x q matrix of elementary modes
%
[m q] = size(EMs);

%%
% Find elementary modes which have a yield higher than the threshold of
% 0.65*maxYield
maxProd = max(input);
thresh = 0.65*maxProd;

%%
% Initialize new EM array
newEMs = [];
iEMs = [];

%%
% Fill new EM array with EMs above threshold
for i=1:q
    if input(i)>= thresh
        newEMs = [newEMs EMs(:,i)];
        iEMs = [iEMs i];
    end
end

%%
% Normalize EMs w.r.t CEL1
normEM = [];

for i=1:size(newEMs,2)
    normEM = newEMs(:,i)/newEMs(1,i);
    newEMs(:,i) = normEM;
end

%%
% Perform SVD to get v_rep of high yield EMs
[U,S,V] = svd(newEMs,'econ');

%%
% Check Principal Component is above SVD threshold
diags = diag(S);
S2 = diags.^2;
S2 = S2/sum(S2);
S2 = cumsum(S2);

% if S2(1)<0.95
%     fprintf('Error: First component of SVD is less than 95% comprehensive.\n')
% elseif S2(2)<0.95
%     fprintf('Error: First two components of SVD are less than 95% comprehensive.\n')
% 
% elseif S2(3)<0.95
%     fprintf('Error: First three components of SVD are less than 95% comprehensive.\n')
%     output = U(:,1:3);
% end

output = U(:,1);
