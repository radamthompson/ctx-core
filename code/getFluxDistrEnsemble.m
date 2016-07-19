function [ ensemble ] = getFluxDistrEnsemble( model, rm_min, rm_max, species, n )
%getFluxDistrEnsemble Make Ensemble of flux distributions
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
%   n := Number of randomly generated flux distributions you wish to
%   calculate
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

% Initialize
rng('shuffle');
ensemble.flux_mat = zeros(length(model.react_name),n);

% Find appropriate dimensions
[model, rand_flux] = getRandomFluxDistr(model, rm_min, rm_max, species );
ensemble.rm = zeros(length(rand_flux.rm),n);
ensemble.calc_rm = zeros(length(rand_flux.calc_rm),n);

% Set initial flux distr
ensemble.flux_mat(:,1) = rand_flux.r;
ensemble.rm(:,1) = rand_flux.rm;
ensemble.calc_rm(:,1) = rand_flux.calc_rm;

% Continue generating ensemble
for i = 2:1:n
    [model, rand_flux] = getRandomFluxDistr(model, rm_min, rm_max, species );
    ensemble.flux_mat(:,i) = rand_flux.r;
    ensemble.rm(:,i) = rand_flux.rm;
    ensemble.calc_rm(:,i) = rand_flux.calc_rm;
end



