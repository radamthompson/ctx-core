function [ angleVec ] = seqFluxAngles( fluxVec1, fluxVec2, react_names )
%seqFluxAngles Sequentially deleting reactions
%   This script compares two flux distributions by sequentially deleting
%   reactions and then finding the angle between the vectors
%
%
% R. Adam Thompson
% Trinh Lab
% University of Tennessee, Knoxville
% rthomp46@utk.edu
%
% Created: May 19, 2014 
% Last edit: 

% Initialize
m1 = length(fluxVec1);
m2 = length(fluxVec2);
angleVec = zeros(m1,1);
angleVecDeg = angleVec;

% Check flux vectors for consistency
if m1~=m2
    error('Flux Distibution Vectors must be the same length.')
end
[r1 c1] = size(fluxVec1);
if r1<= c1
    fluxVec1 = fluxVec1';
end
[r2 c2] = size(fluxVec2);
if r2<= c2
    fluxVec2 = fluxVec2';
end

% Loop to sequentially delete reactions
for i = 1:m1
    if i==1
        newV1 = fluxVec1(2:end);
        newV2 = fluxVec2(2:end);
    end
    if i==m1
        newV1 = fluxVec1(1:end-1);
        newV1 = fluxVec1(1:end-1);
    end
    if i~=1 && i~=m1
        newV1 = [fluxVec1(1:i-1); fluxVec1(i+1:end)];
        newV2 = [fluxVec2(1:i-1); fluxVec2(i+1:end)];
    end
    costheta = dot(newV1,newV2)/(norm(newV1)*norm(newV2));
    angleVec(i) = acos(costheta);
    angleVecDeg(i) = radtodeg(angleVec(i));
end

% Print Table of angles
sortmat = [(1:m1)' angleVec angleVecDeg fluxVec1 fluxVec2];
sortmat = sortrows(sortmat,2);

fprintf('%15s%15s%15s%15s%15s\n','Enzyme','Angle (rad)','Angle (deg)','Flux A','Flux B')
for k = 1:m1
    fprintf('%15s%15.5g%15.5g%15.5g%15.5g\n',react_names{sortmat(k,1)}, sortmat(k,2:end))
end

end
