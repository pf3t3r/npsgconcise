function [ax,p,ks,obs,Sk,Ku,rV,pV,ad,cOut,pOut] = L1_helper(tmp,maxMld,threshold,testSel,hypTest,logAxis,season,suppressFig)
%%L1_helper: this function makes the calculation of KS p-values, skewness,
%%and kurtosis a little more efficient for L1 (the mixed layer). 
% INPUTS
% tmp: the .txt datafile with five columns of data (ID, time, time in
% hours, pressure, and variable concentration). Only the ID, pressure, and
% concentration are used.
% maxMld: the maximum mixed layer depth per cruise. Only values shallower
% than this will be considered.
% unc: uncertainty. To be removed...
% OUTPUTS
% ax = figure identifier, needed for labelling and export,
% p = pressures where sufficient measurements exist in mixed layer,
% ks = KS test p-values at those pressures,
% obs = observations per depth,
% Sk = skewness at depths where ks is taken,
% Ku = kurtosis at the same depths.

if nargin < 8
    suppressFig = false;
end
if nargin < 7
    season = 0;
    % 0 = no seasonal analysis, 1 = winter, 2 = spring, 3 = summer,
    % 4 = autumn
end
if nargin < 6
    logAxis = true;
end
if nargin < 5
    hypTest = "ks";
end
if nargin < 4
    testSel = 4;
end
% Default threshold of 50 based on findings of Mishra et al (2019), Ghasemi
% & Zahediasl (2012), and Ahad et al (2011).
if nargin < 3
    threshold = 50;
end

pIn = tmp.data(:,4);
cIn = tmp.data(:,5);
idIn = tmp.data(:,1);

%%% seasonal analysis
if season ~= 0
    botidIn = tmp.data(:,2);
    n = length(pIn);
    botidIn(botidIn==-9) = nan;
    botId2 = num2str(botidIn);
    botMth = nan(n,1);
    for i = 1:n
        tmpX = str2num(botId2(i,1:end-4));
        if ~isnan(tmpX)
            botMth(i) = tmpX;
        end
    end
    winter = nan(n,1); spring = nan(n,1); summer = nan(n,1); autumn = nan(n,1);
    for i = 1:n
        tmpY = botMth(i);
        if (tmpY == 12) || (tmpY == 1) || (tmpY == 2)
            winter(i) = 1;
        end
        if (tmpY == 3) || (tmpY == 4) || (tmpY == 5)
            spring(i) = 1;
        end
        if (tmpY == 6) || (tmpY == 7) || (tmpY == 8)
            summer(i) = 1;
        end
        if (tmpY == 9) || (tmpY == 10) || (tmpY == 11)
            autumn(i) = 1;
        end
    end

    winIds = []; sprIds = []; sumIds = []; autIds = []; 
    for i = 1:n
        if winter(i) == 1
            winIds = [winIds i];
        end
        if spring(i) == 1
            sprIds = [sprIds i];
        end
        if summer(i) == 1
            sumIds = [sumIds i];
        end
        if autumn(i) == 1
            autIds = [autIds i];
        end
    end

    if season == 1
        cIn = cIn(winIds);
        pIn = pIn(winIds);
        idIn = idIn(winIds);
        %nSea = length(winIds);
    elseif season == 2
        cIn = cIn(sprIds);
        pIn = pIn(sprIds);
        idIn = idIn(sprIds);
        %nSea = length(sprIds);
    elseif season == 3
        cIn = cIn(sumIds);
        pIn = pIn(sumIds);
        idIn = idIn(sumIds);
        %nSea = length(sumIds);
    elseif season == 4
        cIn = cIn(autIds);
        pIn = pIn(autIds);
        idIn = idIn(autIds);
        %nSea = length(autIds);
    end
end
%%% end seasonal analysis



% 2. Extract data in ML
% [idOut,pOut,cOut] = extractMldVals(idIn,pIn,cIn,maxMld);

tmp = num2str(idIn);
bottleCRN = str2num(tmp(:,1:3));
clear tmp;

L = length(pIn);

% OPTION 2
tmpP = nan(1,L);
tmpCrn = nan(1,L);
tmpX = nan(1,L);
tmpId = nan(1,L);

for i = 1:L
    if bottleCRN(i) == 330
        stop = i;
        break
    elseif bottleCRN(i) > 330
        stop = i;
        break
    else
        stop = L+1;
    end
end

% disp(stop);

for i = 1:stop-1
    % OPTION 2: MAX MLD per cruise. This is what we will use.
    tmp = maxMld(bottleCRN(i));
    if pIn(i) < tmp
        tmpP(i) = pIn(i);
        tmpCrn(i) = bottleCRN(i);
        tmpX(i) = cIn(i);
        tmpId(i) = idIn(i);
    end
end

% OPTION 2
pOut = tmpP(~isnan(tmpP));
crnOut = tmpCrn(~isnan(tmpCrn));
cOut = tmpX(~isnan(tmpX));
idOut = tmpId(~isnan(tmpId));


% 3. Bin data
% [~,pOutB,cOutB,~,~] = cleanAndBin(pOut,cOut,idOut');

% Remove bottles that are too close to the surface (< 2.5 dbar)
idRm = pOut > 2.5;
pOut = pOut(idRm);
cOut = cOut(idRm);
botid = idOut(idRm)';

% Remove bottles where concentration of X = 0
idZero = cOut == 0;
pOut = pOut(~idZero);
cOut = cOut(~idZero);
botid = botid(~idZero);

% Save cruise number (CRN) of each bottle - needed below
tmp = num2str(botid);
bottleCRN = str2num(tmp(:,1:3));
clear tmp;

% Remove bottles from cruises 330 on (b/c fluorescence analysis not done)
for i = 1:length(pOut)
    if bottleCRN(i) > 329
        id329 = i - 1;
        break;
    else
        id329 = length(pOut);
    end
end

pOut = pOut(1:id329);
cOut = cOut(1:id329);
clear idRm idZero id329 i;

% pb5 = discretize(pOut,0:5:200);
pb10 = discretize(pOut,0:10:200);
% n5 = max(pb5);
n10 = max(pb10);

cOutB = cOut;
pOutB = pb10;



% 4. Calculate KS p-value, skewness, kurtosis, Vuong Parameters
% [ks,obs,p,Sk,Ku,rV,pV,ad] = ksOfBinnedCon(cOutB,pOutB,10,threshold);

obs = nan(20,1); n = 20; depth = 5:10:200; ad = nan(4,20);
std = nan(1,n);
for i = 1:n
    % find concentration X_i at binned pressure i
    X_i = cOutB(pOutB==i);
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
        ad(:,i) = nan;
        ks(:,i) = nan;
        sk(i) = nan;
        ku(i) = nan;
        rV(:,i) = nan;
        pV(:,i) = nan;
    end
end

tmp = [];
for i = 1:n
    if ~isnan(sum(ad(:,i)))
        tmp = [tmp i];
    end
end
p = depth(tmp);
Sk = sk(tmp);
Ku = ku(tmp);
%sd = sd(tmp,:);
rV = rV(:,tmp);
pV = pV(:,tmp);
ks = ks(:,~all(isnan(ks)));
ad = ad(:,~all(isnan(ad)));


% 4.a. Intercomparison of results from Vuong's Test: easily see best
% distribution at each depth.
vuongRes = nan(1,length(p));

% 4.a.i Default Case.
if testSel == 4
    for i = 1:length(p)
        if rV(1,i) & rV(2,i) & rV(3,i) > 0
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
% 4.a.ii. Normal-Lognormal Case ONLY.
    for i = 1:length(p)
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
if suppressFig == false
    ax = figure;
    plotKs(p,ks,obs,Sk,Ku,0.5,10.5,true,threshold,vuongRes,pV,[0 100],false,hypTest,ad,testSel,"bot",logAxis);
else
    ax = nan;
end

end