function ksVals = getKS(X,n,interval)
%getKS is an intermediate function that call statsplot2 with an optional
% subselection of an input vector.
% INPUT: X = input 2D vector on which Kolmogorov Smirnov test is to be
% applied; n = depth [dbar]; interval = optional subselection, e.g. to only
% look at casts taken in a particular season
% OUTPUT: ksVals = 5 x n vector giving p-values for normal, lognormal,
% weibull, gamma, and exponential distributions. See statsplot2.m for more
% information.

if nargin < 3
    interval = 1:length(X(1,:));
end

for i = 1:n
    tmp = X(i,interval);
    tmp(isnan(tmp)) = [];
    if length(tmp) > 3
        [~,ksVals(:,i),~] = statsplot2(tmp,'noplot');
    end
end
clear tmp;

end