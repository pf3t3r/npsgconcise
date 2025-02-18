function [ax,p,ks,obs,Sk,Ku,rV,pV,ad,cOut,pOut] = L1_helper(tmp,maxMld,threshold,testSel,hypTest,logAxis,season)
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
[idOut,pOut,cOut] = extractMldVals(idIn,pIn,cIn,maxMld);

% 3. Bin data
[~,pOutB,cOutB,~,~] = cleanAndBin(pOut,cOut,idOut');

% 4. Calculate KS p-value, skewness, kurtosis, Vuong Parameters
[ks,obs,p,Sk,Ku,rV,pV,ad] = ksOfBinnedCon(cOutB,pOutB,10,threshold);

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
ax = figure;
plotKs(p,ks,obs,Sk,Ku,0.5,10.5,true,threshold,vuongRes,pV,[0 100],false,hypTest,ad,testSel);

end