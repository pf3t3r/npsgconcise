function [ax,mldCon,rV,pV,ks,vuongRes] = L1_ctdHelper(X,pIn,maxMld,threshold,testSel,hypTest,logAxis)
% [ax,tr,ks,obs,sk,ku,rV,pV]
%%L1_helper: this function makes the calculation of KS p-values, skewness,
%%and kurtosis a little more efficient for L1 (the mixed layer). 
% INPUTS
% tmp: the .txt datafile with five columns of data (ID, time, time in
% hours, pressure, and variable concentration). Only the ID, pressure, and
% concentration are used.
% maxMld: the maximum mixed layer depth per cruise. Only values shallower
% than this will be considered.
% OUTPUTS
% ax = figure identifier, needed for labelling and export,
% p = pressures where sufficient measurements exist in mixed layer,
% ks = KS test p-values at those pressures,
% obs = observations per depth,
% Sk = skewness at depths where ks is taken,
% Ku = kurtosis at the same depths.

% % 'unc' unused for now
% if nargin < 3
%     unc = nan(70,16);
% end

if nargin < 7
    logAxis = true;
end
if nargin < 6
    hypTest = "ks";
end
if nargin < 5
    testSel = 4;
end
% Default threshold = 50 based on findings of Mishra et al (2019), Ghasemi
% & Zahediasl (2012), and Ahad et al (2011).
if nargin < 4
    threshold = 50;
end

% for i = 1:length(epN(1,:))
%     test = maxMld(i);
%     disp(test);
% end

n = length(pIn);

mldCon = nan(size(X));
for k = 1:length(X(1,:))
    for j = 1:n
        if pIn(j) < maxMld(k)
            mldCon(j,k) = X(j,k);
        end
    end
end


% init
ks = nan(5,n); ad = nan(4,n);
rV = nan(10,n); pV = nan(10,n);
sk = nan(1,n); ku = nan(1,n); obs = nan(1,n);

for i = 1:n
    tmp = mldCon(i,:);
    tmp(isnan(tmp) | tmp<0) = 0;
    tmp(tmp==0) = [];
    if length(tmp) > 3
        gammaParams = mle(tmp,"distribution","Gamma");
        pdG = makedist("Gamma",gammaParams(1),gammaParams(2));
        [~,ks(:,i),~] = statsplot2(tmp,'noplot');
        [~,ad(2,i)] = adtest(tmp,"Distribution","logn");
        [~,ad(1,i)] = adtest(tmp,"Distribution","norm");
        [~,ad(3,i)] = adtest(tmp,"Distribution","weibull");
        [~,ad(4,i)] = adtest(tmp,Distribution=pdG,MCTol=0.05);
        [rV(:,i),pV(:,i)] = bbvuong(tmp);
        sk(i) = skewness(tmp);
        ku(i) = kurtosis(tmp);
        obs(i) = length(tmp);
    end
end

% tr = pIn;

for i = 1:n
    if obs(i) < threshold
        ks(:,i) = nan;  
        ad(:,i) = nan;
        sk(i) = nan;
        ku(i) = nan;
        rV(:,i) = nan;
        pV(:,i) = nan;
    end
end

tmp = [];
for i = 1:n
    if ~isnan(sum(ks(:,i)))
        tmp = [tmp i];
    end
end
tr = pIn(tmp);
sk = sk(tmp);
ku = ku(tmp);
rV = rV(:,tmp);
pV = pV(:,tmp);

ks = ks(:,~all(isnan(ks)));
ad = ad(:,~all(isnan(ad)));

% 4.a. Intercomparison of results from Vuong's Test: easily see best
% distribution at each depth.
vuongRes = nan(1,length(tr));

% 4.a.i. 
if testSel == 4 % Default Case.
    for i = 1:length(tr)
        if rV(1,i) > 0 & rV(2,i) > 0 & rV(3,i) > 0
            disp('Normal');
            vuongRes(i) = 1;
        elseif rV(1,i) < 0 & rV(5,i) > 0 & rV(6,i) > 0
            disp('Lognormal');
            vuongRes(i) = 2;
        elseif rV(2,i) < 0 & rV(5,i) < 0 & rV(8,i) > 0
            disp('Weibull');
            vuongRes(i) = 3;
        elseif rV(3,i) < 0 & rV(6,i) < 0 & rV(8,i) < 0
            disp('Gamma');
            vuongRes(i) = 4;
        end
    end
elseif testSel == 2
    % Normal vs Lognormal Case.
    for i = 1:length(tr)
        if rV(1,i) > 0
            disp('Normal');
            vuongRes(i) = 1;
        elseif rV(1,i) < 0
            disp('Lognormal');
            vuongRes(i) = 2;
        end
    end
end

% 5. Plot results
ax = figure;
plotKs(tr,ks,obs,sk,ku,0,100,true,threshold,vuongRes,pV,[0 100],true,hypTest,ad,testSel,"ctd",logAxis);
% 
% disp(vuongRes);
end