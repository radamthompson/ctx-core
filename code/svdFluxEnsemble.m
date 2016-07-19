function [ ensemble ] = svdFluxEnsemble( ensemble )
%svdFluxEnsemble Perform SVD on flux ensemble
%   This function will take a flux ensemble structure from
%   getFluxDistrEnsemble() and perform the SVD decomposition of the flux
%   ensemble.
%
%   model := Metatool model structure, requires elementary modes to be in
%   field model.ems
%   ensemble := Flux ensemble structure obtained from running
%   getFluxDistrEnsemble()
%
%   Output:
%   ensemble := Ensemble structure with SVD stats
%
% R. Adam Thompson
% Trinh Lab
% University of Tennessee, Knoxville
% rthomp46@utk.edu
%
% Created: October 6, 2015 
% Last edit: 

vs = real(ensemble.flux_mat);
vssums = sqrt(sum(vs.^2));
vs = vs*diag(1./vssums);

% Perform SVD
[ensemble.Uv ensemble.Sv ensemble.Vv] = svd(vs,'econ');

% Calculate mean to compare to first principal component vector from U
ensemble.vmean = mean(vs,2);
% Calculate standard deviation for each of the fluxes
ensemble.vstd = std(vs,0,2);

% Construct scree plot
Sdiag = diag(ensemble.Sv);
S2 = Sdiag.^2;
S2 = S2/sum(S2);
S2 = cumsum(S2);
figure
plot(S2)
xlabel('Number of principal components')
ylabel('Fraction of cumulative sum of squared singular values')
grid on
xlim([0 10])

ensemble.S2 = S2;

