function [ model ] = runMetatool( fname )
%runMetatool Take a Metatool input file and calculate elementary modes
%   fname := Text file with metabolic network set up as described in
%   Metatool documentation
%
%   output := structure of Metatool network and elementary modes
%
% R. Adam Thompson
% Trinh Lab
% University of Tennessee, Knoxville
% rthomp46@utk.edu
%
% Created: February 12, 2012
% Last edit: February 10, 2015


% Bring in model
ex = parse(fname);

% Send to metatool
model = metatool(ex);

% Elementary modes for non-reduced network
model.ems = model.sub'*model.rd_ems;

% Stochiometric matrix
model.stoich = model.ext*model.ems;

end

