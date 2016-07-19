function [ theta, newVec ] = projectEMs( ems, vector )
%projectEMs 
%   This script takes in a set of elementary modes and projects a flux
%   vector onto the space spanned by said modes. Then, it calculates the
%   angle in degrees between the projection and the input flux vector

%   Written by: R. Adam Thompson, UTK, 5/19/2014
%   Last Modified: 5/23/14 

ProjMat = ems*ems';

newVec = ProjMat*vector;

costheta = dot(newVec,vector) / (norm(newVec) * norm(vector));
theta = acos(costheta);
theta = radtodeg(theta);

end

