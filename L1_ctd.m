% Script to output L1 ctd results for the statistical analysis.

close all; clc; clear;
addpath("baroneRoutines\"); addpath("func\"); addpath("output\");
set(groot, "defaultFigureUnits", "centimeters", "defaultFigurePosition", [3 3 15 15]);

% Test Cases
principalAnalysisAd = true;     % main analysis with A-D test
principleAnalysisKs = false;    % main analysis with K-S test
seasonalAnalysis = false;       % seasonality of statistics
noOfDists = 2;
adThresh = 30;                  % only for chl-a in A-D
logAxes = true;                 % output p-values as log values (true)
if logAxes == true
    lp = "log/";
else
    lp = "";
end

%% Load MaxMld and Chl-a (EpN) and CTD data

% F2 = 131:329
epN = load("output\CTD\chla.mat").meanEpN;
epN = epN(1:76,131:329);
pIn = 0:2:150;
maxMld = load("mldVals.mat").maxMld;

ctdData = load("datafiles\ctd_iso_ALL.mat").ctd;
msng = [21, 48, 207, 218, 276];
cR = 1:1:329;
cRm = setdiff(cR,msng);

% Temperature T
T = nan(329,76,31);
meanT = nan(76,329);
for i = cRm
    tmp = ctdData(i).t(1:76,:);
    if length(tmp) > 3
        for j = 1:length(tmp(1,:))
            T(i,:,j) = tmp(:,j);
        end
    end
end
for i = 1:329
    meanT(:,i) = mean(squeeze(T(i,:,:)),2,"omitnan");
end

% Salinity SP
% S missing data same as for T? Yes!
SP = nan(329,76,31);
meanSp = nan(76,329);
for i = cRm
    tmp = ctdData(i).sp(1:76,:);
    if length(tmp) > 3
        for j = 1:length(tmp(1,:))
            SP(i,:,j) = tmp(:,j);
        end
    end
end
for i = 1:329
    meanSp(:,i) = mean(squeeze(SP(i,:,:)),2,"omitnan");
end

% O2
O2 = nan(329,76,31);
meanO2 = nan(76,329);
for i = cRm
    tmp = ctdData(i).o(1:76,:);
    if length(tmp) > 3
        for j = 1:length(tmp(1,:))
            O2(i,:,j) = tmp(:,j);
        end
    end
end
for i = 1:329
    meanO2(:,i) = mean(squeeze(O2(i,:,:)),2,"omitnan");
end

%% GSW Processing.
% Convert temperature and practical salinity into conservative temperature
% and absolute salinity, and then calculate potential density.

stnALOHA_lon = -158;
stnALOHA_lat = 22.75;
SA = gsw_SA_from_SP(meanSp,pIn',stnALOHA_lon,stnALOHA_lat);
CT = gsw_CT_from_t(SA,meanT,pIn');
sigma0 = gsw_sigma0(SA,CT);
[N2,p_mid] = gsw_Nsquared(SA,CT,pIn',stnALOHA_lat);

%% Show Derived Parameters

% Potential Density Anomaly (wrt p=0) 'sigma0'
figure;
scatter(sigma0,pIn,5,"yellow");
set(gca,"YDir","reverse");

% Absolute Salinity SA
figure;
scatter(SA,pIn,4,"yellow");
set(gca,"YDir","reverse");

% Conservative Temperature CT
figure;
scatter(CT,pIn,4,"yellow");
set(gca,"YDir","reverse");

%% Show vertical profiles

% figure
% plot(meanT,pIn); set(gca,"YDir","reverse");
% title("Temperature (C): 88-21");
% 
% figure
% plot(meanSp,pIn); set(gca,"YDir","reverse");
% title("Practical Salinity (g/kg): 88-21");
% 
% figure
% plot(meanO2,pIn); set(gca,"YDir","reverse");
% title("$O_2$ (mmol m$^{-3}$)",Interpreter="latex");

%% Principal Analysis
if principleAnalysisKs == true
    %  K-S
    tmpT = "";
    
    % CHL-A
    ax = L1_ctdHelper(epN,pIn,maxMld,50,noOfDists);
    sgtitle("fluorescence (01-21): L1");
    exportgraphics(ax,"figures/L1/ctd/"+lp+"chla" + tmpT + ".png"); clear ax;
    % save("output\L1\ctd\chla.mat","p","ks","obs","Sk","Ku");
    
    % T
    [ax,mldCon,rV,pV,ks,vuongRes] = L1_ctdHelper(meanT,pIn,maxMld,50,noOfDists);
    sgtitle("CTD Temperature 88-21: L1");
    exportgraphics(ax,"figures/L1/ctd/"+lp+"T" + tmpT + ".png"); clear ax;
    % save("output\L1\ctd\t.mat","p","ks","obs","Sk","Ku");
    
    % SP
    ax = L1_ctdHelper(meanSp,pIn,maxMld,50,noOfDists);
    sgtitle("CTD S_P 88-21: L1");
    exportgraphics(ax,"figures/L1/ctd/"+lp+"sp" + tmpT + ".png"); clear ax;
    % save("output\L1\ctd\sp.mat","p","ks","obs","Sk","Ku");
    
    % O2
    [ax,mldCon] = L1_ctdHelper(meanO2,pIn,maxMld,50,noOfDists);
    sgtitle("CTD O2 88-21: L1");
    exportgraphics(ax,"figures/L1/ctd/"+lp+"o2" + tmpT + ".png"); clear ax;
    % save("output\L1\ctd\o2.mat","p","ks","obs","Sk","Ku");
    
end

if principalAnalysisAd == true

    % A-D
    tmpT = "-ad";
    
    tmpx = "";
    % CHL-A
    ax = L1_ctdHelper(epN,pIn,maxMld,adThresh,noOfDists,"ad");
    sgtitle("L1"+tmpx);
    exportgraphics(ax,"figures/L1/ctd/"+lp+"chla" + tmpT + ".png"); clear ax;
    % save("output\L1\ctd\chla.mat","p","ks","obs","Sk","Ku");
    
    % T
    ax = L1_ctdHelper(meanT,pIn,maxMld,adThresh,noOfDists,"ad");
    sgtitle("CTD Temperature 88-21: L1"+tmpT);
    exportgraphics(ax,"figures/L1/ctd/"+lp+"T" + tmpT + ".png"); clear ax;
    % save("output\L1\ctd\t.mat","p","ks","obs","Sk","Ku");

    % CT
    ax = L1_ctdHelper(CT,pIn,maxMld,adThresh,noOfDists,"ad");
    sgtitle("L1");
    exportgraphics(ax,"figures/L1/ctd/"+lp+"CT" + tmpT + ".png"); clear ax;
    % save("output\L1\ctd\t.mat","p","ks","obs","Sk","Ku");
    
    % SP
    [ax,mldCon,rV,pV,ad,V] = L1_ctdHelper(meanSp,pIn,maxMld,adThresh,noOfDists,"ad");
    sgtitle("CTD S_P 88-21: L1"+tmpT);
    exportgraphics(ax,"figures/L1/ctd/"+lp+"sp" + tmpT + ".png"); clear ax;
    % save("output\L1\ctd\sp.mat","p","ks","obs","Sk","Ku");

    % SA
    ax = L1_ctdHelper(SA,pIn,maxMld,adThresh,noOfDists,"ad");
    sgtitle("CTD S_A 88-21: L1"+tmpT);
    exportgraphics(ax,"figures/L1/ctd/"+lp+"sa" + tmpT + ".png"); clear ax;
    % save("output\L1\ctd\sp.mat","p","ks","obs","Sk","Ku");

    % sigma0
    ax = L1_ctdHelper(sigma0,pIn,maxMld,adThresh,noOfDists,"ad");
    sgtitle("CTD \sigma_0 88-21: L1"+tmpT);
    exportgraphics(ax,"figures/L1/ctd/"+lp+"sigma0" + tmpT + ".png"); clear ax;
    % save("output\L1\ctd\sp.mat","p","ks","obs","Sk","Ku");

    % N^2 Stratification
    ax = L1_ctdHelper(N2,p_mid(:,1),maxMld,adThresh,noOfDists,"ad");
    sgtitle("N^2 88-21: L1"+tmpT);
    exportgraphics(ax,"figures/L1/ctd/"+lp+"n_squared" + tmpT + ".png"); clear ax;
    % save("output\L1\ctd\sp.mat","p","ks","obs","Sk","Ku");
    
    
    % O2
    ax = L1_ctdHelper(meanO2,pIn,maxMld,adThresh,noOfDists,"ad");
    sgtitle("CTD O2 88-21: L1"+tmpT);
    exportgraphics(ax,"figures/L1/ctd/"+lp+"o2" + tmpT + ".png"); clear ax;
    % save("output\L1\ctd\o2.mat","p","ks","obs","Sk","Ku");
end

%% Set-up seasonal test.

% find "average" month of a cruise, then give an ID to that cruise saying
% which season it is (Spring, Summer, Autumn, or Winter).

avgMonth = nan(329,1);
winter = nan(329,1);
spring = nan(329,1);
summer = nan(329,1);
autumn = nan(329,1);

for i = 1:329
    tmp = round(mean(month(datetime(ctdData(i).date,"ConvertFrom","datenum"))));
    
    if (tmp == 12) || (tmp == 1) || (tmp == 2)
        winter(i) = 1;
    end
    if (tmp == 3) || (tmp == 4) || (tmp == 5)
        spring(i) = 1;
    end
    if (tmp == 6) || (tmp == 7) || (tmp == 8)
        summer(i) = 1;
    end
    if (tmp == 9) || (tmp == 10) || (tmp == 11)
        autumn(i) = 1;
    end

    avgMonth(i) = tmp;
end

% ASTRO seasons -> winter = jfm; spring = amj; summer = jas; autumn = ond.

winIds = []; sprIds = []; sumIds = []; autIds = []; 
for i = 1:329
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

% Calculate average mld per season
mldW = mean(maxMld(winIds),"omitnan");
mldSp = mean(maxMld(sprIds),"omitnan");
mldSu = mean(maxMld(sumIds),"omitnan");
mldA = mean(maxMld(autIds),"omitnan");

% Average dcm per season
dcm = load("dcm.mat").meanPcm;  
dcmW = mean(dcm(winIds),"omitnan");
dcmSp = mean(dcm(sprIds),"omitnan");
dcmSu = mean(dcm(sumIds),"omitnan");
dcmA = mean(dcm(autIds),"omitnan");

%% Test Seasonal Analysis

% Naming scheme
% -01 = winter, -02 = spring, -03 = summer, -04 = autumn

if seasonalAnalysis == true

    % K-S
    % winter
    tmpT = "-01";

    chla = load("output\CTD\chla.mat").meanEpN(1:101,winIds(36:end));
    ax = L1_ctdHelper(chla,pIn,maxMld,50,noOfDists,"ks");
    sgtitle("chl-a " + tmpT);
    exportgraphics(ax,"figures/L1/ctd/"+lp+"chla" + tmpT + ".png");

    % T
    ax = L1_ctdHelper(meanT(:,winIds),pIn,maxMld,50,noOfDists,"ks");
    sgtitle("T " + tmpT);
    exportgraphics(ax,"figures/L1/ctd/"+lp+"T" + tmpT + ".png");
    
    % Sp
    ax = L1_ctdHelper(meanSp(:,winIds),pIn,maxMld,50,noOfDists,"ks");
    sgtitle("$S_p$ " + tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L1/ctd/" + lp + "Sp" + tmpT + ".png");
    
    % O2
    ax = L1_ctdHelper(meanO2(:,winIds),pIn,maxMld,50,noOfDists,"ks");
    sgtitle("$O_2$ " + tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L1/ctd/" + lp + "O2" + tmpT + ".png");

    % spring
    tmpT = "-02";

    chla = load("output\CTD\chla.mat").meanEpN(1:101,sprIds(33:end));
    ax = L1_ctdHelper(chla,pIn,maxMld,50,noOfDists,"ks");
    sgtitle("chl-a " + tmpT);
    exportgraphics(ax,"figures/L1/ctd/"+lp+"chla" + tmpT + ".png");

    % T
    ax = L1_ctdHelper(meanT(:,sprIds),pIn,maxMld,50,noOfDists,"ks");
    sgtitle("T " + tmpT);
    exportgraphics(ax,"figures/L1/ctd/"+lp+"T" + tmpT + ".png");
    
    % Sp
    ax = L1_ctdHelper(meanSp(:,sprIds),pIn,maxMld,50,noOfDists,"ks");
    sgtitle("$S_p$ " + tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L1/ctd/" + lp + "Sp" + tmpT + ".png");
    
    % O2
    ax = L1_ctdHelper(meanO2(:,sprIds),pIn,maxMld,50,noOfDists,"ks");
    sgtitle("$O_2$ " + tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L1/ctd/" + lp + "O2" + tmpT + ".png");

    % summer
    tmpT = "-03";

    chla = load("output\CTD\chla.mat").meanEpN(1:101,sumIds(33:end));
    ax = L1_ctdHelper(chla,pIn,maxMld,50,noOfDists,"ks");
    sgtitle("chl-a " + tmpT);
    exportgraphics(ax,"figures/L1/ctd/"+lp+"chla" + tmpT + ".png");

    % T
    ax = L1_ctdHelper(meanT(:,sumIds),pIn,maxMld,50,noOfDists,"ks");
    sgtitle("T " + tmpT);
    exportgraphics(ax,"figures/L1/ctd/"+lp+"T" + tmpT + ".png");
    
    % Sp
    ax = L1_ctdHelper(meanSp(:,sumIds),pIn,maxMld,50,noOfDists,"ks");
    sgtitle("$S_p$ " + tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L1/ctd/" + lp + "Sp" + tmpT + ".png");
    
    % O2
    ax = L1_ctdHelper(meanO2(:,sumIds),pIn,maxMld,50,noOfDists,"ks");
    sgtitle("$O_2$ " + tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L1/ctd/" + lp + "O2" + tmpT + ".png");

    % autumn
    tmpT = "-04";

    chla = load("output\CTD\chla.mat").meanEpN(1:101,sumIds(30:end));
    ax = L1_ctdHelper(chla,pIn,maxMld,50,noOfDists,"ks");
    sgtitle("chl-a " + tmpT);
    exportgraphics(ax,"figures/L1/ctd/"+lp+"chla" + tmpT + ".png");

    % T
    ax = L1_ctdHelper(meanT(:,autIds),pIn,maxMld,50,noOfDists,"ks");
    sgtitle("T " + tmpT);
    exportgraphics(ax,"figures/L1/ctd/"+lp+"T" + tmpT + ".png");
    
    % Sp
    ax = L1_ctdHelper(meanSp(:,autIds),pIn,maxMld,50,noOfDists,"ks");
    sgtitle("$S_p$ " + tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L1/ctd/" + lp + "Sp" + tmpT + ".png");
    
    % O2
    ax = L1_ctdHelper(meanO2(:,autIds),pIn,maxMld,50,noOfDists,"ks");
    sgtitle("$O_2$ " + tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L1/ctd/" + lp + "O2" + tmpT + ".png");

    % A-D
    % winter
    tmpT = "-ad-01";

    chla = load("output\CTD\chla.mat").meanEpN(1:101,winIds(36:end));
    ax = L1_ctdHelper(chla,pIn,maxMld,adThresh,noOfDists,"ad");
    sgtitle("chl-a " + tmpT);
    exportgraphics(ax,"figures/L1/ctd/"+lp+"chla" + tmpT + ".png");

    % T
    ax = L1_ctdHelper(meanT(:,winIds),pIn,maxMld,adThresh,noOfDists,"ad");
    sgtitle("T " + tmpT);
    exportgraphics(ax,"figures/L1/ctd/"+lp+"T" + tmpT + ".png");
    
    % Sp
    ax = L1_ctdHelper(meanSp(:,winIds),pIn,maxMld,adThresh,noOfDists,"ad");
    sgtitle("$S_p$ " + tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L1/ctd/" + lp + "Sp" + tmpT + ".png");
    
    % O2
    ax = L1_ctdHelper(meanO2(:,winIds),pIn,maxMld,adThresh,noOfDists,"ad");
    sgtitle("$O_2$ " + tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L1/ctd/" + lp + "O2" + tmpT + ".png");

    % spring
    tmpT = "-ad-02";

    chla = load("output\CTD\chla.mat").meanEpN(1:101,sprIds(33:end));
    ax = L1_ctdHelper(chla,pIn,maxMld,adThresh,noOfDists,"ad");
    sgtitle("chl-a " + tmpT);
    exportgraphics(ax,"figures/L1/ctd/"+lp+"chla" + tmpT + ".png");

    % T
    ax = L1_ctdHelper(meanT(:,sprIds),pIn,maxMld,adThresh,noOfDists,"ad");
    sgtitle("T " + tmpT);
    exportgraphics(ax,"figures/L1/ctd/"+lp+"T" + tmpT + ".png");
    
    % Sp
    ax = L1_ctdHelper(meanSp(:,sprIds),pIn,maxMld,adThresh,noOfDists,"ad");
    sgtitle("$S_p$ " + tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L1/ctd/" + lp + "Sp" + tmpT + ".png");
    
    % O2
    ax = L1_ctdHelper(meanO2(:,sprIds),pIn,maxMld,adThresh,noOfDists,"ad");
    sgtitle("$O_2$ " + tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L1/ctd/" + lp + "O2" + tmpT + ".png");

    % summer
    tmpT = "-ad-03";

    chla = load("output\CTD\chla.mat").meanEpN(1:101,sumIds(33:end));
    ax = L1_ctdHelper(chla,pIn,maxMld,adThresh,noOfDists,"ad");
    sgtitle("chl-a " + tmpT);
    exportgraphics(ax,"figures/L1/ctd/"+lp+"chla" + tmpT + ".png");

    % T
    ax = L1_ctdHelper(meanT(:,sumIds),pIn,maxMld,adThresh,noOfDists,"ad");
    sgtitle("T " + tmpT);
    exportgraphics(ax,"figures/L1/ctd/"+lp+"T" + tmpT + ".png");
    
    % Sp
    ax = L1_ctdHelper(meanSp(:,sumIds),pIn,maxMld,adThresh,noOfDists,"ad");
    sgtitle("$S_p$ " + tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L1/ctd/" + lp + "Sp" + tmpT + ".png");
    
    % O2
    ax = L1_ctdHelper(meanO2(:,sumIds),pIn,maxMld,adThresh,noOfDists,"ad");
    sgtitle("$O_2$ " + tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L1/ctd/" + lp + "O2" + tmpT + ".png");

    % autumn
    tmpT = "-ad-04";

    chla = load("output\CTD\chla.mat").meanEpN(1:101,sumIds(30:end));
    ax = L1_ctdHelper(chla,pIn,maxMld,adThresh,noOfDists,"ad");
    sgtitle("chl-a " + tmpT);
    exportgraphics(ax,"figures/L1/ctd/"+lp+"chla" + tmpT + ".png");

    % T
    ax = L1_ctdHelper(meanT(:,autIds),pIn,maxMld,adThresh,noOfDists,"ad");
    sgtitle("T " + tmpT);
    exportgraphics(ax,"figures/L1/ctd/"+lp+"T" + tmpT + ".png");
    
    % Sp
    ax = L1_ctdHelper(meanSp(:,autIds),pIn,maxMld,adThresh,noOfDists,"ad");
    sgtitle("$S_p$ " + tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L1/ctd/" + lp + "Sp" + tmpT + ".png");
    
    % O2
    ax = L1_ctdHelper(meanO2(:,autIds),pIn,maxMld,adThresh,noOfDists,"ad");
    sgtitle("$O_2$ " + tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L1/ctd/" + lp + "O2" + tmpT + ".png");
end


%% Unused.
% % Load Nitrate (NO3-)
% 
% NO3 = nan(329,76,31);
% meanNO3 = nan(76,329);
% 
% for i = cRm
%     tmp = ctdData(i).n(1:76,:);
%     if length(tmp) > 3
%         for j = 1:length(tmp(1,:))
%             NO3(i,:,j) = tmp(:,j);
%         end
%     end
% end
% 
% for i = 1:329
%     meanNO3(:,i) = mean(squeeze(NO3(i,:,:)),2,"omitnan");
% end
% %
% figure
% plot(meanNO3,pIn); set(gca,"YDir","reverse");
% title("$NO_3^{-}$ (mmol m$^{-3}$)",Interpreter="latex");
% %% NO3-
% [ax,p,ks,obs,Sk,Ku,rV,pV] = L1_ctdHelper(meanNO3,pIn,maxMld);
% sgtitle("CTD NO3- 88-21: L1");
% exportgraphics(ax,"figures/L1/ctd/notUsed/no3" + tmpT + ".png"); clear ax;
% % save("output\L1\ctd\no3.mat","p","ks","obs","Sk","Ku");