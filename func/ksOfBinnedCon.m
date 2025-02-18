function [ks,obs,tr,sk,ku,rV,pV,ad] = ksOfBinnedCon(X, p, binning, threshold)
%ksOfBinnedCon find the KS statistic
% INPUTS:
% X = substance concentration,
% p = binned pressure [dbar]
% binning = range of depths to bin values (default=10)
% threshold = no. of values needed for us to consider results (default=100)
% OUTPUTS: 
% ks = KS test for five distributions 
% obs = observations per depth
% tr = array of depths above threshold (=100)
% sk = skewness
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
    ks = nan(5,40); obs = nan(40,1); n = 40; depth = 5:5:200; ad = nan(4,40);
elseif binning == 10
    ks = nan(5,20); obs = nan(20,1); n = 20; depth = 5:10:200; ad = nan(4,20);
else
    msg = 'Binning input not valid. Must be either 5 or 10 dbar.';
    error(msg);
end

std = nan(1,n);
for i = 1:n
    % find concentration X_i at binned pressure i
    X_i = X(p==i);
    % apply KS test to X_i
    % change limit below to >3 to fix error with picoeu -> may change other
    % results    
    if length(X_i) > 3
        gammaParams = mle(X_i,"distribution","Gamma");
        pdG = makedist("Gamma",gammaParams(1),gammaParams(2));
        [~,ks(:,i),~] = statsplot2(X_i,'noplot');
        [~,ad(2,i)] = adtest(X_i,"Distribution","logn","Alpha",0.005);
        [~,ad(1,i)] = adtest(X_i,"Distribution","norm","Alpha",0.005);
        [~,ad(3,i)] = adtest(X_i,"Distribution","weibull");
        [~,ad(4,i)] = adtest(X_i,Distribution=pdG,MCTol=0.05);
        [rV(:,i),pV(:,i)] = bbvuong(X_i);
        sk(i) = skewness(X_i);
        ku(i) = kurtosis(X_i);
    end
    obs(i) = length(X_i);
    clear X_i;
end

for i = 1:n
    if obs(i) < threshold
        ks(:,i) = nan;
        ad(:,i) = nan;
        sk(i) = nan;
        ku(i) = nan;
        %sd(i,:) = nan;
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
tr = depth(tmp);
sk = sk(tmp);
ku = ku(tmp);
%sd = sd(tmp,:);
rV = rV(:,tmp);
pV = pV(:,tmp);
ks = ks(:,~all(isnan(ks)));
ad = ad(:,~all(isnan(ad)));

end