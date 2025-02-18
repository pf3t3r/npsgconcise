function [trange,ks,obsPerBin,Sk,Ku,bottleArray,sigB,Xout] = ksOfIsoLagrangian(id,sig,dcmArray,X,threshold)
%ksOfIsoLagrangian quickly find the DCM-centred (Lagrangian) transformation for a
%given variable in binned isopycnal (density) coordinates.
% INPUTS:
% id: bottle ID
% sig = binned densities
% dcmArray = shows where the DCM is
% X = bottle concentration
% Ltid = the 'mystery variable'
% OUTPUTS:
% bottleArray = not needed anymore, to remove
% trange = depths where sufficient measurements are present
% Xout = also not used??
% sigB = binned density, also not used?
% ks = Kolmogorov-Smirnov Statistic (p value). High p-value indicates that
% the given distribution fits the data better.
% obsPerBin = no. of observations in a particular depth bin. 


if nargin < 5
    threshold = 100;
end

CRN = str2num(id(:,1:3)); 
cast = str2num(id(:,6:8));
cast(cast==100) = nan;

bottleArray = [CRN cast sig];

% Save unique cruise/cast combinations
t = rmmissing(unique(bottleArray(:,1:2),"rows"));

% Find array location of DCM at cruise & cast required.
dcmArrayRowNo = [];
for i = 1:length(dcmArray(:,1))
    for x = 1:length(t)
        if dcmArray(i,1:2) == t(x,1:2) 
            dcmArrayRowNo = [dcmArrayRowNo i];
        end
    end
end

% Split bottle concentration by cruise & cast
tid = [];
for i = 2:length(sig)
    if bottleArray(i,1) > bottleArray(i-1,1) || bottleArray(i,2) > bottleArray(i-1,2)
        tid = [tid i];
    end
end

tSigCm = nan(length(sig),1);
tSigCm(1:tid(1)-1) = dcmArray(dcmArrayRowNo(1),4);
tSigCm(tid(end):end) = dcmArray(dcmArrayRowNo(end),4);

% The number of unique cruise & cast combinations
Ltid = length(tid);
disp(length(tid));
for i = 2:Ltid-2
    %disp(i);
    tSigCm(tid(i):tid(i+1)-1) = dcmArray(dcmArrayRowNo(i),4);
end

tSigCmR = round(tSigCm,3,"significant");
bottleArray = [bottleArray tSigCm tSigCmR];

tPLagrangian = bottleArray(:,3) - bottleArray(:,5);
bottleArray = [bottleArray tPLagrangian];

bottleArray = [bottleArray X];
tmin = min(bottleArray(:,6));
tmax = max(bottleArray(:,6));
trange = tmin:0.1:tmax;

trange = round(trange,3);

Xout = bottleArray(:,7);
sigB = bottleArray(:,6);
ks = nan(5,length(trange));
obsPerBin = nan(1,length(trange));

Sk = nan(1,length(trange)); Ku = nan(1,length(trange));

for i = 1:length(trange)
    % Load bottle measurements at each binned density level
    tmp2 = round(10000*sigB);
    tmp = Xout(tmp2 == round(10000*trange(i)));
    tmp(tmp<=0) = nan;
    tmp(isnan(tmp)) = [];
    obsPerBin(i) = length(tmp);
    if length(tmp) > 3
        disp(i);
        [~,ks(:,i),~,~,~,~] = statsplot2(tmp,'noplot');
        Sk(i) = skewness(tmp);
        Ku(i) = kurtosis(tmp);
    end
end

for i = 1:length(trange)
    if obsPerBin(i) < threshold
        ks(:,i) = nan;
        Sk(i) = nan;
        Ku(i) = nan;
    end
end

end