function Check = EmExPa(ems,st,irrev_react)
%EmExPa Check if EMs are ExPas
% This function takes a set of elementary modes and determines which are
% extreme pathways. 
%
%
% R. Adam Thompson
% Trinh Lab
% University of Tennessee, Knoxville
% rthomp46@utk.edu
%

n_ems = size(ems,2);
irrev_ndx = find(irrev_react == 1);
n_irrev = length(irrev_ndx);
rev_ndx = find(irrev_react == 0);
n_rev = length(rev_ndx);


for i = 1:n_ems
    zbar = find(ems(:,i)~=0);
    zbarR = union(zbar,rev_ndx);
    EM_check = length(zbar)-rank(st(:,zbar));
    ExPa_check = length(zbarR)-rank(st(:,zbarR));
    if EM_check == 1
        Check.EM(i,1) = 1;
    else
        Check.EM(i,1) = 0;
    end
    if ExPa_check == 1;
        Check.ExPa(i,1) = 1;
    else
        Check.ExPa(i,1) = 0;
    end
    Check.ExPa_check(i,1) = ExPa_check;
end
