function [tr,ks,obs,sk,ku,rV,pV,ad,X_out,bottleArray,pB] = ksOfLagrangian(id,p,dcm,X,threshold)
%ksOfLagrangian(): quickly find the DCM-centred (Lagrangian) transformation for a
%given variable.
% INPUTS:
% id: bottle ID
% p = pressure
% dcm = shows where the DCM is
% X = bottle concentration
% Ltid = not sure how this works but hey it works...
% OUTPUTS:
% bottleArray = not needed anymore, to remove
% tr = depths where sufficient measurements are present
% X_out = also not used??
% pB = binned pressure, also not used?
% ks = Kolmogorov-Smirnov Statistic (p value). High p-value indicates that
% the given distribution fits the data better.
% obs = no. of observations in a particular depth bin. 

if nargin < 5
    threshold = 50;
end

% Extract cruise number 'crn' and 'cast'
crn = str2num(id(:,1:3)); 
cast = str2num(id(:,6:8));
cast(cast==100) = nan;

bottleArray = [crn cast p];

% Create an array of all unique bottle cruise/cast combinations
botCrnCast = rmmissing(unique(bottleArray(:,1:2),"rows"));

% Find these unique cruise/cast combinations in the 'dcm' array
dcmCrnCast = [];
for i = 1:length(dcm(:,1))
    for x = 1:length(botCrnCast)
        if dcm(i,1:2) == botCrnCast(x,1:2) 
            dcmCrnCast = [dcmCrnCast i];
        end
    end
end

% Split bottle concentration by cruise & cast
tid = [];
for i = 2:length(p)
    % check if CRN or CAST changes
    if bottleArray(i,1) > bottleArray(i-1,1) || bottleArray(i,2) > bottleArray(i-1,2)
        tid = [tid i];
    end
end

tPcm = nan(length(p),1);
tPcm(1:tid(1)-1) = dcm(dcmCrnCast(1),3);
tPcm(tid(end):end) = dcm(dcmCrnCast(end),3);

Ltid = length(tid);

tmp = length(dcmCrnCast) - 1;
Ltid = tmp;
for i = 2:Ltid-2
    tPcm(tid(i):tid(i+1)-1) = dcm(dcmCrnCast(i),3);
end

bottleArray = [bottleArray tPcm];

% THIS is where we convert to Lagrangian pressure coordinates!!!
tPLagrangian = nan(length(p),1);
tPLagrangian = bottleArray(:,3) - bottleArray(:,4);
bottleArray = [bottleArray tPLagrangian];

pB10 = round(bottleArray(:,5),-1);
bottleArray = [bottleArray pB10];

bottleArray = [bottleArray X];
tmin = min(bottleArray(:,6));
tmax = max(bottleArray(:,6));
tr = tmin:10:tmax;

X_out = bottleArray(:,7);
pB = bottleArray(:,6);
ks = nan(5,length(tr));
ad = nan(4,length(tr));
rV = nan(10,length(tr));
pV = nan(10,length(tr));
obs = nan(1,length(tr));

sk = nan(1,length(tr));

for i = 1:length(tr)
    tmp = X_out(pB==tr(i));
    tmp(tmp<=0) = nan;
    tmp(isnan(tmp)) = [];
    obs(i) = length(tmp);
    if length(tmp) > 3
        gammaParams = mle(tmp,"distribution","Gamma");
        pdG = makedist("Gamma",gammaParams(1),gammaParams(2));
        [~,ks(:,i),~] = statsplot2(tmp,'noplot');
        [~,ad(2,i)] = adtest(tmp,"Distribution","logn","Alpha",0.005);
        [~,ad(1,i)] = adtest(tmp,"Distribution","norm","Alpha",0.005);
        [~,ad(3,i)] = adtest(tmp,"Distribution","weibull");
        [~,ad(4,i)] = adtest(tmp,Distribution=pdG,MCTol=0.05);
        [rV(:,i),pV(:,i)] = bbvuong(tmp);
        sk(i) = skewness(tmp);
        ku(i) = kurtosis(tmp);
%         tmpDat = [std(tmp) std(log(tmp))];
%         tmpDatMu = [mean(tmp) mean(log(tmp))];
%         tmpComp = [tmpMle(1)/tmpDat(1) tmpMle(2)/tmpDat(2)];
%         tmpCompMu = [muMle(1)/tmpDatMu(1) muMle(2)/tmpDatMu(2)];
%         mu(i,:) = [muMle(1) muMle(2) tmpDatMu(1) tmpDatMu(2) tmpCompMu(1) tmpCompMu(2)];
%         sd(i,:) = [tmpMle(1) tmpMle(2) tmpDat(1) tmpDat(2) tmpComp(1) tmpComp(2)];
%         c95(i,1:2) = [tmpC95(1,1) tmpC95(2,1)];
%         c95(i,3:4) = [tmpC95(1,2) tmpC95(2,2)];
%         c95(i,5:6) = [tmpC95(1,3) tmpC95(2,3)];
%         c95(i,7:8) = [tmpC95(1,4) tmpC95(2,4)];
    end
end

for i = 1:length(tr)
    if obs(i) < threshold
        ks(:,i) = nan;
        ad(:,i) = nan;
        rV(:,i) = nan;
        sk(i) = nan;
        ku(i) = nan;
%         sd(i,:) = nan;
%         c95(i,:) = nan;
%         mu(i,:) = nan;
    end
end

% tmp = [];
% for i = 1:length(tr)
%     if ~isnan(sum(ks(:,i)))
%         tmp = [tmp i];
%     end
% end
% 
% trange2 = tr(tmp);
% sk = sk(tmp);
% ku = ku(tmp);
% sd = sd(tmp,:);
% ks = ks(:,~all(isnan(ks)));

end