function [ax,X_out,pSubml,bA,ks,obs,sk,ku,rV,p,pV,ad,pr,vuongRes,cSubml,pB] = L2_helper(tmp,maxMld,dcm,threshold,testSel,hypTest,yLimits,yLimitsObs,season)
%%L2_helper: this function makes the calculation of KS p-values, skewness,
%%and kurtosis a little more efficient for L2 (sub-mixed layer region that
% is centred on the DCM). 
% INPUTS
% tmp: the .txt datafile with five columns of data (ID, time, time in
% hours, pressure, and variable concentration). Only the ID, pressure, and
% concentration are used.
% maxMld: the maximum mixed layer depth per cruise. Only values shallower
% than this will be considered.
% dcm: pressure of deep chlorophyll maximum (by cruise)
% OUTPUTS
% ax = figure identifier, needed for labelling and export,
% p = pressures where sufficient measurements exist in mixed layer,
% ks = KS test p-values at those pressures,
% obs = observations per depth,
% Sk = skewness at depths where ks is taken,
% Ku = kurtosis at the same depths.

logAxis = true; % true => output p-values in log x-axis, otherwise no log plot.

% 1. Assign default parameter values.
% Default threshold of 50 based on Mishra et al (2019), Ghasemi & Zahediasl 
% (2012), and Ahad et al (2011).
if nargin < 4
    threshold = 50;
end
% Default to test against four distributions in K-S or A-D.
if nargin < 5
    testSel = 4;
end
% Default = "ks" (alternative = "ad").
if nargin < 6
    hypTest = "ks";
end
% Default y-limits for K-S/A-D and Vuong plots.
if nargin < 7
    yLimits = [-60 60];
end
% Default y-limits for horizontal bar chart showing no. of observations per
% depth. This should vary in a similar way to yLimits above. Test run the
% code to output "pr" to check this.
if nargin < 8
    yLimitsObs = [7 19];
end
if nargin < 9
    season = 0;
    % 0 = no seasonal analysis, 1 = winter, 2 = spring, 3 = summer,
    % 4 = autumn
end

id = num2str(tmp.data(:,1));
p = tmp.data(:,4);
c = tmp.data(:,5);
% clear tmp;

idOld = id;
%%% seasonal analysis
if season ~= 0
    botidIn = tmp.data(:,2);
    n = length(p);
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
        c = c(winIds);
        p = p(winIds);
        id = id(winIds,:);
        %nSea = length(winIds);
    elseif season == 2
        c = c(sprIds);
        p = p(sprIds);
        id = id(sprIds,:);
        %nSea = length(sprIds);
    elseif season == 3
        c = c(sumIds);
        p = p(sumIds);
        id = id(sumIds,:);
        %nSea = length(sumIds);
    elseif season == 4
        c = c(autIds);
        p = p(autIds);
        id = id(autIds,:);
        %nSea = length(autIds);
    end
end
%%% end seasonal analysis

% 2. Extract data beneath ML
[idSubml,pSubml,cSubml] = extractSMLC(id,p,c,maxMld);

% 3. Calculate KS p-value, skewness, kurtosis
% ..., centre around DCM (?)
[pr,ks,obs,sk,ku,rV,pV,ad,X_out,bA,pB] = ksOfLagrangian(idSubml,pSubml,dcm,cSubml,threshold);

% 3.a. Intercomparison of results from Vuong's Test: easily see best
% distribution at each depth.
vuongRes = zeros(1,length(pr));
rV(isnan(rV)) = 0;

if testSel==4
    for i = 1:length(pr)
        %disp(i);
        if rV(1,i) & rV(2,i) & rV(3,i) > 0
            %disp('Normal');
            vuongRes(i) = 1;
        elseif rV(1,i) < 0 & rV(5,i) > 0 & rV(6,i) > 0
            %disp('Lognormal');
            vuongRes(i) = 2;
        elseif rV(2,i) < 0 & rV(5,i) < 0 & rV(8,i) > 0
            %disp('Weibull');
            vuongRes(i) = 3;
        elseif rV(3,i) < 0 & rV(6,i) < 0 & rV(8,i) < 0
            %disp('Gamma');
            vuongRes(i) = 4;
        end
    end
    rV(rV==0) = nan;
elseif testSel == 2
    for i = 1:length(pr)
        if rV(1,i)  > 0
            vuongRes(i) = 1;
        elseif rV(1,i) < 0
            vuongRes(i) = 2;
        end
    end
    rV(rV==0) = nan;
end

% limits = [pr(barchartLimits(1)) pr(barchartLimits(2))];
limits = yLimits;
obsId = [yLimitsObs(1) yLimitsObs(2)];

% 4. Plot results
ax = figure;
[a,b] = plotKs2(pr,ks,obs,sk,ku,limits,threshold,vuongRes,obsId,pV,hypTest,ad,testSel,logAxis);

% disp(vuongRes);
end