function [pKs, pLil, pAd, pSw, obs, depth2, Sk, ku] = pValueOfBinned(X, p, binning, threshold)
%pValueOfBinned finds the p-value corresponding to a concentration
%following three statistical tests: the Kolmogorov-Smirnov test (K-S), the
%Lilliefors-corrected K-S test (Lil), and the Anderson-Darling test (A-D).
% INPUTS:
% X = substance concentration,
% p = binned pressure [dbar]
% binning = range of depths to bin values (default=10)
% threshold = no. of values needed for us to consider results (default=50)
% OUTPUTS: 
% pKs = p-value (K-S)
% pLil = p-value (Lil)
% pAd = p-value (A-D)
% obs = observations per depth
% depth2 = array of depths above threshold (=100)
% Sk = skewness
% ku = kurtosis
% sd = standard deviation (3 x n). 1: STD of mle (norm); 2: STD of data
% (norm); 3: stdMle/stdData (norm)

if nargin <4
    threshold = 50;
end

if nargin <3
    binning = 10;
end

if binning == 5
    pKs = nan(5,40); pLil = nan(5,40); pAd = nan(5,40); pSw = nan(2,40); obs = nan(40,1); n = 40; depth = 5:5:200;
elseif binning == 10
    pKs = nan(5,20); pLil = nan(5,20); pAd = nan(5,20); pSw = nan(2,20); obs = nan(20,1); n = 20; depth = 5:10:200;
else
    msg = 'Binning input not valid. Must be either 5 or 10 dbar.';
    error(msg);
end

for i = 1:n
    % find concentration X_i at binned pressure i
    X_i = X(p==i);
    % apply KS test to X_i
    % change limit below to >3 to fix error with picoeu -> may change other
    % results
    if length(X_i) > 3
        [~,pKs(:,i),pLil(:,i),pAd(:,i),pSw(:,i)] = statsplotComp(X_i);
        %disp(size(tmpC95));
        Sk(i) = skewness(X_i);
        ku(i) = kurtosis(X_i);
    end
    obs(i) = length(X_i);
    clear X_i;
end

for i = 1:n
    if obs(i) < threshold
        pKs(:,i) = nan; pLil(:,i) = nan; pAd(:,i) = nan; pSw(:,i) = nan;
        Sk(i) = nan;
        ku(i) = nan;
    end
end

tmp = [];
% tmp2 = [];
for i = 1:n
    if ~isnan(sum(pKs(:,i)))
        tmp = [tmp i];
    end
%     if ~isnan(sum(pLil(:,i)))
%         tmp2 = [tmp2 i];
%     end
%     if ~isnan(sum(pAd(:,i)))
%         tmp3 = [tmp3 i];
%     end
end
depth2 = depth(tmp);
Sk = Sk(tmp);
ku = ku(tmp);

pKs = pKs(:,~all(isnan(pKs)));
pLil = pLil(:,~all(isnan(pLil)));
pAd = pAd(:,~all(isnan(pAd)));
pSw = pSw(:,~all(isnan(pSw)));

end