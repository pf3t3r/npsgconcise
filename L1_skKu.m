% Skewness-kurtosis analysis of the mixed layer (L1) at Station ALOHA for chl-a, other pigments, and BGC variables.
% We import the data, process it, and calculate the skewness and kurtosis.

% This script duplicates a lot of L1_bot right now. Must clean this up.

clear; clc; close all;
addpath("baroneRoutines\"); addpath("func\");
set(groot, "defaultFigureUnits", "centimeters", "defaultFigurePosition", [3 3 10 15]);


%% Variables
% Set these parameters each time the script is run.

season = 0;     % if zero, run default analysis. otherwise run seasonal
                % analysis (1 = winter, 2 = spring, 3 = summer, 4 =
                % autumn)
threshold = 30; % default threshold for statistical tests

%% Rest of code.

%% Import MLD.
ctdData = importdata("datafiles\ctd_iso_ALL.mat").ctd;
maxMld = nan(329,1);
for i = 1:329
    if ~isnan([ctdData(i).mld003])
        maxMld(i) = max([ctdData(i).mld003]);
    end
end

clear ctdData i;
save mldVals.mat maxMld;

%% Import Data
tmp = importdata("data/L1/hplcChla_88-21_150.txt");
pIn = tmp.data(:,4);
cIn = tmp.data(:,5);
idIn = tmp.data(:,1);

%% seasonal analysis code
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

%% Extract data in mixed layer.
% 2. Extract data in ML
% [idOut,pOut,cOut] = extractMldVals(idIn,pIn,cIn,maxMld);
% Extract the cruise number 'crn'.
tmp = num2str(idIn);
crn = str2num(tmp(:,1:3)); clear tmp;

L = length(pIn);    % Length of dataset

% Initialise ...
tmpP = nan(1,L);
tmpCrn = nan(1,L);
tmpX = nan(1,L);
tmpId = nan(1,L);

% Stop at cruise number 329. This was in order to match the latest data
% available for CTD fluorometry but I could consider updating this.
for i = 1:L
    if crn(i) == 330
        stop = i;
        break
    elseif crn(i) > 330
        stop = i;
        break
    else
        stop = L+1;
    end
end

% Save values at pressures above (i.e. shallower than) the mixed layer depth
for i = 1:stop-1
    tmp = maxMld(crn(i));
    if pIn(i) < tmp
        tmpP(i) = pIn(i);
        tmpCrn(i) = crn(i);
        tmpX(i) = cIn(i);
        tmpId(i) = idIn(i);
    end
end

% Remove NaNs from new arrays of values in mixed layer.
pOut = tmpP(~isnan(tmpP));
crnOut = tmpCrn(~isnan(tmpCrn));
cOut = tmpX(~isnan(tmpX));
idOut = tmpId(~isnan(tmpId));

%% Data: Quality Control and Bin.

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
crn = str2num(tmp(:,1:3)); clear tmp;

% Remove bottles from cruises 330 on (b/c fluorescence analysis not done)
for i = 1:length(pOut)
    if crn(i) > 329
        id329 = i - 1;
        break;
    else
        id329 = length(pOut);
    end
end

pOut = pOut(1:id329);
cOut = cOut(1:id329);
clear idRm idZero id329 i;

pb10 = discretize(pOut,0:10:200);
n10 = max(pb10);

cOutB = cOut;
pOutB = pb10;

%% Calculate skewness and kurtosis.
% Calculate the skewness and kurtosis of the cleaned data. NOTE that other
% parameters related to A-D calculation are still present. These are to be
% removed.

obs = nan(20,1); n = 20; depth = 5:10:200; ad = nan(4,20);
for i = 1:n
    % find concentration X_i at binned pressure i
    X_i = cOutB(pOutB==i);
    % apply statistical tests to the data   
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

% Remove values that do not meet the threshold
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

% Remove these nan values
tmp = [];
for i = 1:n
    if ~isnan(sum(ad(:,i)))
        tmp = [tmp i];
    end
end
p = depth(tmp);
Sk = sk(tmp);
Ku = ku(tmp);
rV = rV(:,tmp);
pV = pV(:,tmp);
ks = ks(:,~all(isnan(ks)));
ad = ad(:,~all(isnan(ad)));

% 4.a. Intercomparison of results from Vuong's Test: easily see best
% distribution at each depth.
vuongRes = nan(1,length(p));

testSel = 2;
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

n = length(p);

% % Create Annotations for Vuong's Test Results
% annot = strings(1,n);
% anClr = strings(1,n);
% anClr(cellfun(@isempty,anClr)) = '#FFFFFF';
% tmpEmph = strings(1,n); tmpEmph(cellfun(@isempty,tmpEmph)) = 'bold';
% % Default Case
% 
% alphaHy = 0.005;
% alphaLlr = 0.10;

% if testSel == 4
%     for i = 1:n
%         if strcmp(hypTest,"ks")
%             if vuongRes(i) == 1 && ks(1,i) > alphaHy
%                 % Remove label if only one dist is not rejected by K-S.
%                 if length(find(ks(:,i)>alphaHy)) == 1
%                     tmp = "";
%                 else
%                     tmp = "Normal";
%                 end
%                 anClr(i) = '#a6cee3';
%                 if pV(1,i) > alphaLlr && ks(2,i) > alphaHy
%                     tmpEmph(i) = 'normal';
%                     tmp = append(tmp," L");
%                 end
%                 if pV(2,i) > alphaLlr && ks(3,i) > alphaHy
%                     tmpEmph(i) = 'normal';
%                     tmp = append(tmp," W");
%                 end
%                 if pV(3,i) > alphaLlr && ks(4,i) > alphaHy
%                     tmpEmph(i) = 'normal';
%                     tmp = append(tmp," G");
%                 end
%                 annot(i) = tmp;
%             elseif vuongRes(i) == 2 && ks(2,i) > alphaHy
%                 % Remove label if only one dist is not rejected by K-S.
%                 if length(find(ks(:,i)>alphaHy)) == 1
%                     tmp = "";
%                 else
%                     tmp = "Lognormal";
%                 end
%                 anClr(i) = '#1f78b4';
%                 if pV(1,i) > alphaLlr && ks(1,i) > alphaHy
%                     tmpEmph(i) = 'normal';
%                     tmp = append(tmp," N");
%                 end
%                 if pV(5,i) > alphaLlr && ks(3,i) > alphaHy
%                     tmpEmph(i) = 'normal';
%                     tmp = append(tmp," W");
%                 end
%                 if pV(6,i) > alphaLlr && ks(4,i) > alphaHy
%                     tmpEmph(i) = 'normal';
%                     tmp = append(tmp," G");
%                 end
%                 annot(i) = tmp;
%             elseif vuongRes(i) == 3 && ks(3,i) > alphaHy
%                 % Remove label if only one dist is not rejected by K-S.
%                 if length(find(ks(:,i)>alphaHy)) == 1
%                     tmp = "";
%                 else
%                     tmp = "Weibull";
%                 end
%                 anClr(i) = '#b2df8a';
%                 if pV(2,i) > alphaLlr && ks(1,i) > alphaHy
%                     tmpEmph(i) = 'normal';
%                     tmp = append(tmp," N");
%                 end
%                 if pV(5,i) > alphaLlr && ks(2,i) > alphaHy
%                     tmpEmph(i) = 'normal';
%                     tmp = append(tmp," L");
%                 end
%                 if pV(8,i) > alphaLlr && ks(4,i) > alphaHy
%                     tmpEmph(i) = 'normal';
%                     tmp = append(tmp," G");
%                 end
%                 annot(i) = tmp;
%             elseif vuongRes(i) == 4 && ks(4,i) > alphaHy
%                 % Remove label if only one dist is not rejected by K-S.
%                 if length(find(ks(:,i)>alphaHy)) == 1
%                     tmp = "";
%                 else
%                     tmp = "Gamma";
%                 end
%                 anClr(i) = '#33a02c';
%                 if pV(6,i) > alphaLlr && ks(2,i) > alphaHy
%                     tmpEmph(i) = 'normal';
%                     tmp = append(tmp," L");
%                 end
%                 if pV(3,i) > alphaLlr && ks(1,i) > alphaHy
%                     tmpEmph(i) = 'normal';
%                     tmp = append(tmp," N");
%                 end
%                 if pV(8,i) > alphaLlr && ks(3,i) > alphaHy
%                     tmpEmph(i) = 'normal';
%                     tmp = append(tmp," W");
%                 end
%                 annot(i) = tmp;
%             elseif vuongRes(i) == 0
%                 annot(i) = "";
%             end
%         elseif strcmp(hypTest,"ad")
%             % A-D
%             if vuongRes(i) == 1 && ad(1,i) > alphaHy
%                 % Remove label if only one dist is not rejected by K-S.
%                 if length(find(ad(:,i)>alphaHy)) == 1
%                     tmp = "";
%                 else
%                     tmp = "Normal";
%                 end
%                 anClr(i) = '#a6cee3';
%                 if pV(1,i) > alphaLlr && ad(2,i) > alphaHy
%                     tmpEmph(i) = 'normal';
%                     tmp = append(tmp," L");
%                 end
%                 if pV(2,i) > alphaLlr && ad(3,i) > alphaHy
%                     tmpEmph(i) = 'normal';
%                     tmp = append(tmp," W");
%                 end
%                 if pV(3,i) > alphaLlr && ad(4,i) > alphaHy
%                     tmpEmph(i) = 'normal';
%                     tmp = append(tmp," G");
%                 end
%                 annot(i) = tmp;
%             elseif vuongRes(i) == 2 && ad(2,i) > alphaHy
%                 % Remove label if only one dist is not rejected by K-S.
%                 if length(find(ad(:,i)>alphaHy)) == 1
%                     tmp = "";
%                 else
%                     tmp = "Lognormal";
%                 end
%                 anClr(i) = '#1f78b4';
%                 if pV(1,i) > alphaLlr && ad(1,i) > alphaHy
%                     tmpEmph(i) = 'normal';
%                     tmp = append(tmp," N");
%                 end
%                 if pV(5,i) > alphaLlr && ad(3,i) > alphaHy
%                     tmpEmph(i) = 'normal';
%                     tmp = append(tmp," W");
%                 end
%                 if pV(6,i) > alphaLlr && ad(4,i) > alphaHy
%                     tmpEmph(i) = 'normal';
%                     tmp = append(tmp," G");
%                 end
%                 annot(i) = tmp;
%             elseif vuongRes(i) == 3 && ad(3,i) > alphaHy
%                 % Remove label if only one dist is not rejected by K-S.
%                 if length(find(ad(:,i)>alphaHy)) == 1
%                     tmp = "";
%                 else
%                     tmp = "Weibull";
%                 end
%                 anClr(i) = '#b2df8a';
%                 if pV(2,i) > alphaLlr && ad(1,i) > alphaHy
%                     tmpEmph(i) = 'normal';
%                     tmp = append(tmp," N");
%                 end
%                 if pV(5,i) > alphaLlr && ad(2,i) > alphaHy
%                     tmpEmph(i) = 'normal';
%                     tmp = append(tmp," L");
%                 end
%                 if pV(8,i) > alphaLlr && ad(4,i) > alphaHy
%                     tmpEmph(i) = 'normal';
%                     tmp = append(tmp," G");
%                 end
%                 annot(i) = tmp;
%             elseif vuongRes(i) == 4 && ad(4,i) > alphaHy
%                 % Remove label if only one dist is not rejected by K-S.
%                 if length(find(ad(:,i)>alphaHy)) == 1
%                     tmp = "";
%                 else
%                     tmp = "Gamma";
%                 end
%                 anClr(i) = '#33a02c';
%                 if pV(6,i) > alphaLlr && ad(2,i) > alphaHy
%                     tmpEmph(i) = 'normal';
%                     tmp = append(tmp," L");
%                 end
%                 if pV(3,i) > alphaLlr && ad(1,i) > alphaHy
%                     tmpEmph(i) = 'normal';
%                     tmp = append(tmp," N");
%                 end
%                 if pV(8,i) > alphaLlr && ad(3,i) > alphaHy
%                     tmpEmph(i) = 'normal';
%                     tmp = append(tmp," W");
%                 end
%                 annot(i) = tmp;
%             elseif vuongRes(i) == 0
%                 annot(i) = "";
%             end
% 
%         end
%     end
% elseif testSel == 2
%     % Normal-Lognormal Case
%     for i = 1:n
%         if strcmp(hypTest,"ks")
%             if vuongRes(i) == 1 && ks(1,i) > alphaHy
%                 annot(i) = "Normal";
%                 anClr(i) = '#c51b7d';
%                 if pV(1,i) > alphaLlr
%                     tmpEmph(i) = 'normal';
%                 end
%             elseif vuongRes(i) == 2 && ks(2,i) > alphaHy
%                 annot(i) = "Lognormal";
%                 anClr(i) = '#4d9221';
%                 if pV(1,i) > alphaLlr
%                     tmpEmph(i) = 'normal';
%                 end
%             elseif vuongRes(i) == 0
%                 annot(i) = "";
%             end
%         elseif strcmp(hypTest,"ad")
%             if vuongRes(i) == 1 && ad(1,i) > alphaHy
%                 annot(i) = "Normal";
%                 anClr(i) = '#c51b7d';
%                 if pV(1,i) > alphaLlr
%                     tmpEmph(i) = 'normal';
%                 end
%             elseif vuongRes(i) == 2 && ad(2,i) > alphaHy
%                 annot(i) = "Lognormal";
%                 anClr(i) = '#4d9221';
%                 if pV(1,i) > alphaLlr
%                     tmpEmph(i) = 'normal';
%                 end
%             elseif vuongRes(i) == 0
%                 annot(i) = "";
%             end
%         end
%     end
% end

% Lognormal family: generate theoretical skewness and kurtosis
sigTh = linspace(0,1,1000);
for i = 1:length(sigTh)
    skLogn(i) = (exp(sigTh(i)^2) + 2)*(sqrt(exp(sigTh(i)^2) - 1));
    kuLogn(i) = exp(4*sigTh(i)^2) + 2*exp(3*sigTh(i)^2) + 3*exp(2*sigTh(i)^2) - 3;
end

% Negative Distributions
skLognN = -skLogn;
kuLognN = kuLogn;


kurtLimB = 10; skewLimA = 0; skewLimB = 2.5;
if max(ku) > 10 & min(sk) < 0
    kurtLimB = max(ku) + 1;
    skewLimA = min(sk) - 0.1;
    skewLimB = max(sk) + 0.1;
elseif max(ku) > 10
    kurtLimB = max(ku) + 1;
    skewLimB = max(sk) + 0.1;
elseif min(sk) < 0 
    skewLimA = min(sk) - 0.1;
elseif max(sk) > 2.5
    kurtLimB = max(ku) + 1;
    skewLimB = max(sk) + 0.1;
else 
    kurtLimB = 10;
    skewLimA = 0;
    skewLimB = 2.5;
end

scatter(nan,nan,72,[0.8 0.8 0.8],DisplayName='Data');
hold on
scatter(0,3,72,[0.2 0.2 0.2],'DisplayName','Normal',Marker='pentagram',LineWidth=2.5);
if testSel == 2
    plot(skLogn,kuLogn,'DisplayName','Lognormal','Color','#808080',LineStyle='-',LineWidth=1.3);
    plot(skLognN,kuLognN,'Color','#808080',LineStyle='-',LineWidth=1.3,HandleVisibility='off');
elseif testSel == 4
    plot(skLogn,kuLogn,'DisplayName','Logn.','Color','#1f78b4',LineStyle='--',LineWidth=1.7);
    plot(skLognN,kuLognN,'Color','#1f78b4',LineStyle='--',LineWidth=1.7,HandleVisibility='off');
    plot(skWbl,kuWbl,'DisplayName','Weib.','Color','#b2df8a',LineStyle='-',LineWidth=1.7);
    plot(skWblN,kuWblN,'Color','#b2df8a',LineStyle='-',LineWidth=1.7,HandleVisibility='off');
    plot(skGam,kuGam,'DisplayName','Gam.','Color','#33a02c',LineStyle='--',LineWidth=1.7);
    plot(skGamN,kuGamN,'Color','#33a02c',LineStyle='--',LineWidth=1.7,HandleVisibility='off');
    scatter(2,9,'DisplayName','Exp.',Marker='+',LineWidth=1);
    scatter(0,9/5,'DisplayName','Uni.',Marker='*',LineWidth=1);
    scatter(0,21/5,'DisplayName','Logi.',Marker='.',LineWidth=1);
    scatter(1.1395,5.4,'DisplayName','LEV',Marker='x',LineWidth=1);
end
scatter(Sk,Ku,72,[0.8 0.8 0.8],HandleVisibility="off");
clr = 1:1:length(p);
scatter(Sk,Ku,54,clr,"filled","o",HandleVisibility="off");
colormap(gca,cbrewer2("Greens"));
cbar = colorbar;
cbar.Direction = "reverse";
cbar.Ticks = 1:1:length(p);
% cbar.TickLabels = p(1):10:p(end);
cbar.TickLabels = p;
cbar.Label.String = "P [dbar]";
cbar.Label.Position = [0.7 1-0.35];
cbar.Label.Rotation = 0;
% hold on
% add polynomial
% [skS,id] = sort(sk);
% kuS = ku(id);
% [p,S] = polyfit(skS,kuS,2);
% [f,delta] = polyval(p,skS,S);
% plot(skS,f,'r-',DisplayName="Fit");
% plot(skS,f+2*delta,'m--',DisplayName='95% Prediction Interval');
% plot(skS,f-2*delta,'m--',HandleVisibility='off');
hold off
grid minor;
ylim([1 kurtLimB]); xlim([skewLimA skewLimB]);
xlabel('Skewness','FontSize',13,'Interpreter','latex'); ylabel('Kurtosis',FontSize=13,Interpreter='latex');
lgd = legend('Location','best');
% title(lgd,'Distributions');
title('L1','Interpreter','latex','FontSize',13);
% sgtitle("L2 chl-$a$ skewness-kurtosis 1988-2021","Interpreter","latex");
exportgraphics(ax,"figures/L1/bottle/log/chla_ad_skKu.png");