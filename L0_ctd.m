% Script to output L0 ctd results for the statistical analysis.

clear; clc; close all;
addpath("baroneRoutines\");
addpath("func\");
set(groot, 'defaultFigureUnits', 'centimeters', 'defaultFigurePosition', [3 3 28 15]);

% Test Cases
principleAnalysis = true;       % main analysis
startYearAnalysis = false;      % effect of altering start year on HTs
seasonalAnalysis = false;       % seasonality of statistics
logAxes = true;                 % output p-values as log values (true)
if logAxes == true
    lp = "log/";
else
    lp = "";
end

%% Load CTD data: T, Sp, O2.
% CTD data is mean of all casts. Downloaded via FTP and averaged in other
% file. It SHOULD be equivalent to the CTD chloropigment we find on the
% HOT-DOGS system.

ctdData = load("datafiles\ctd_iso_ALL.mat").ctd;
msng = [21, 48, 207, 218, 276];
cR = 1:1:329;
cRm = setdiff(cR,msng);

% Temperature T
T = nan(329,101,31);
meanT = nan(101,329);
for i = cRm
    tmp = ctdData(i).t(1:101,:);
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
SP = nan(329,101,31);
meanSp = nan(101,329);
for i = cRm
    tmp = ctdData(i).sp(1:101,:);
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
O2 = nan(329,101,31);
meanO2 = nan(101,329);
for i = cRm
    tmp = ctdData(i).o(1:101,:);
    if length(tmp) > 3
        for j = 1:length(tmp(1,:))
            O2(i,:,j) = tmp(:,j);
        end
    end
end
for i = 1:329
    meanO2(:,i) = mean(squeeze(O2(i,:,:)),2,"omitnan");
end

%% Derived Quantities

tmpP = 0:2:200;
lon = -158; lat = 22.75;

% Absolute Salinity SA
SA = gsw_SA_from_SP(meanSp,tmpP',lon,lat);

% Conservative Temperature CT
CT = gsw_CT_from_t(SA,meanT,tmpP');

%% Check Parameters for CHLA

chla = load("output\CTD\chla.mat").meanEpN(1:101,131:329);
chla = 1000*chla; % convert to ng/l

muStar = nan(1,101);
sigStar = nan(1,101);
for i = 1:101
    tmp = mle(chla(i,:),'distribution','Lognormal');
    muStar(i) = tmp(1);
    sigStar(i) = tmp(2);
end
depths = 0:2:200;

sigStar = exp(sigStar);
muStar = exp(muStar);

figure
subplot(1,2,1)
plot(muStar,depths); set(gca,"YDir","reverse"); ylabel("P (dbar)");
title("$\mu^*$",Interpreter="latex");
subplot(1,2,2)
plot(sigStar,depths); set(gca,"YDir","reverse");
title("$\sigma^*$",Interpreter="latex");
sgtitle("chlorophyll fluorescence parameters: L0 (2001-2021)");

% sigStar seems to low. Test parameters of random lognormal distributions.
test = lognrnd(-0.6361,0.2654,1000,1);
phatT = mle(test,distribution="Lognormal");
figure
histogram(test);

[h,p,adstat,cv] = adtest(test,"Distribution","logn");

%% Principal Analysis

if principleAnalysis == true

    % K-S
    tmpT = "";
    
    % chla
    chla = load("output\CTD\chla.mat").meanEpN(1:101,131:329);
    ax = L0_ctdHelper(chla,"ks",logAxes);
    sgtitle("chl-a " + tmpT);
    exportgraphics(ax,"figures/L0/ctd/"+lp+"chla" + tmpT + ".png");
    
    % T
    ax = L0_ctdHelper(meanT,"ks",logAxes);
    sgtitle("T " + tmpT);
    exportgraphics(ax,"figures/L0/ctd/"+lp+"T" + tmpT + ".png");
    
    % Sp
    ax = L0_ctdHelper(meanSp,"ks",logAxes);
    sgtitle("$S_p$ " + tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/ctd/" + lp + "Sp" + tmpT + ".png");
    
    % O2
    ax = L0_ctdHelper(meanO2,"ks",logAxes);
    sgtitle("$O_2$ " + tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/ctd/" + lp + "O2" + tmpT + ".png");
    
    % A-D
    tmpT = "-ad";
    
    % chla
    chla = load("output\CTD\chla.mat").meanEpN(1:101,131:329);
    ax = L0_ctdHelper(chla,"ad",logAxes);
    sgtitle("chl-a " + tmpT);
    exportgraphics(ax,"figures/L0/ctd/" + lp + "chla" + tmpT + ".png");
    
    % T
    ax = L0_ctdHelper(meanT,"ad",logAxes);
    sgtitle("T " + tmpT);
    exportgraphics(ax,"figures/L0/ctd/" + lp + "T" + tmpT + ".png");

    % CT
    ax = L0_ctdHelper(CT,"ad",logAxes);
    sgtitle("CT ");
    exportgraphics(ax,"figures/L0/ctd/" + lp + "CT" + tmpT + ".png");

    % Sp
    ax = L0_ctdHelper(meanSp,"ad",logAxes);
    sgtitle("$S_p$ " + tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/ctd/" + lp + "Sp" + tmpT + ".png");
    
    % O2
    ax = L0_ctdHelper(meanO2,"ad",logAxes);
    sgtitle("$O_2$ " + tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/ctd/" + lp + "O2" + tmpT + ".png");

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

%% Test Seasonal Analysis

% Naming scheme
% -01 = winter, -02 = spring, -03 = summer, -04 = autumn

if seasonalAnalysis == true

    % K-S
    % winter
    tmpT = "-01";

    chla = load("output\CTD\chla.mat").meanEpN(1:101,winIds(36:end));
    ax = L0_ctdHelper(chla,"ks",logAxes);
    sgtitle("chl-a " + tmpT);
    exportgraphics(ax,"figures/L0/ctd/"+lp+"chla" + tmpT + ".png");

    % T
    ax = L0_ctdHelper(meanT(:,winIds),"ks",logAxes);
    sgtitle("T " + tmpT);
    exportgraphics(ax,"figures/L0/ctd/"+lp+"T" + tmpT + ".png");
    
    % Sp
    ax = L0_ctdHelper(meanSp(:,winIds),"ks",logAxes);
    sgtitle("$S_p$ " + tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/ctd/" + lp + "Sp" + tmpT + ".png");
    
    % O2
    ax = L0_ctdHelper(meanO2(:,winIds),"ks",logAxes);
    sgtitle("$O_2$ " + tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/ctd/" + lp + "O2" + tmpT + ".png");

    % spring
    tmpT = "-02";

    chla = load("output\CTD\chla.mat").meanEpN(1:101,sprIds(33:end));
    ax = L0_ctdHelper(chla,"ks",logAxes);
    sgtitle("chl-a " + tmpT);
    exportgraphics(ax,"figures/L0/ctd/"+lp+"chla" + tmpT + ".png");

    % T
    ax = L0_ctdHelper(meanT(:,sprIds),"ks",logAxes);
    sgtitle("T " + tmpT);
    exportgraphics(ax,"figures/L0/ctd/"+lp+"T" + tmpT + ".png");
    
    % Sp
    ax = L0_ctdHelper(meanSp(:,sprIds),"ks",logAxes);
    sgtitle("$S_p$ " + tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/ctd/" + lp + "Sp" + tmpT + ".png");
    
    % O2
    ax = L0_ctdHelper(meanO2(:,sprIds),"ks",logAxes);
    sgtitle("$O_2$ " + tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/ctd/" + lp + "O2" + tmpT + ".png");

    % summer
    tmpT = "-03";

    chla = load("output\CTD\chla.mat").meanEpN(1:101,sumIds(33:end));
    ax = L0_ctdHelper(chla,"ks",logAxes);
    sgtitle("chl-a " + tmpT);
    exportgraphics(ax,"figures/L0/ctd/"+lp+"chla" + tmpT + ".png");

    % T
    ax = L0_ctdHelper(meanT(:,sumIds),"ks",logAxes);
    sgtitle("T " + tmpT);
    exportgraphics(ax,"figures/L0/ctd/"+lp+"T" + tmpT + ".png");
    
    % Sp
    ax = L0_ctdHelper(meanSp(:,sumIds),"ks",logAxes);
    sgtitle("$S_p$ " + tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/ctd/" + lp + "Sp" + tmpT + ".png");
    
    % O2
    ax = L0_ctdHelper(meanO2(:,sumIds),"ks",logAxes);
    sgtitle("$O_2$ " + tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/ctd/" + lp + "O2" + tmpT + ".png");

    % autumn
    tmpT = "-04";

    chla = load("output\CTD\chla.mat").meanEpN(1:101,sumIds(30:end));
    ax = L0_ctdHelper(chla,"ks",logAxes);
    sgtitle("chl-a " + tmpT);
    exportgraphics(ax,"figures/L0/ctd/"+lp+"chla" + tmpT + ".png");

    % T
    ax = L0_ctdHelper(meanT(:,autIds),"ks",logAxes);
    sgtitle("T " + tmpT);
    exportgraphics(ax,"figures/L0/ctd/"+lp+"T" + tmpT + ".png");
    
    % Sp
    ax = L0_ctdHelper(meanSp(:,autIds),"ks",logAxes);
    sgtitle("$S_p$ " + tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/ctd/" + lp + "Sp" + tmpT + ".png");
    
    % O2
    ax = L0_ctdHelper(meanO2(:,autIds),"ks",logAxes);
    sgtitle("$O_2$ " + tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/ctd/" + lp + "O2" + tmpT + ".png");

    % A-D
    % winter
    tmpT = "-ad-01";

    chla = load("output\CTD\chla.mat").meanEpN(1:101,winIds(36:end));
    ax = L0_ctdHelper(chla,"ad",logAxes);
    sgtitle("chl-a " + tmpT);
    exportgraphics(ax,"figures/L0/ctd/"+lp+"chla" + tmpT + ".png");

    % T
    ax = L0_ctdHelper(meanT(:,winIds),"ad",logAxes);
    sgtitle("T " + tmpT);
    exportgraphics(ax,"figures/L0/ctd/"+lp+"T" + tmpT + ".png");
    
    % Sp
    ax = L0_ctdHelper(meanSp(:,winIds),"ad",logAxes);
    sgtitle("$S_p$ " + tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/ctd/" + lp + "Sp" + tmpT + ".png");
    
    % O2
    ax = L0_ctdHelper(meanO2(:,winIds),"ad",logAxes);
    sgtitle("$O_2$ " + tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/ctd/" + lp + "O2" + tmpT + ".png");

    % spring
    tmpT = "-ad-02";

    chla = load("output\CTD\chla.mat").meanEpN(1:101,sprIds(33:end));
    ax = L0_ctdHelper(chla,"ad",logAxes);
    sgtitle("chl-a " + tmpT);
    exportgraphics(ax,"figures/L0/ctd/"+lp+"chla" + tmpT + ".png");

    % T
    ax = L0_ctdHelper(meanT(:,sprIds),"ad",logAxes);
    sgtitle("T " + tmpT);
    exportgraphics(ax,"figures/L0/ctd/"+lp+"T" + tmpT + ".png");
    
    % Sp
    ax = L0_ctdHelper(meanSp(:,sprIds),"ad",logAxes);
    sgtitle("$S_p$ " + tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/ctd/" + lp + "Sp" + tmpT + ".png");
    
    % O2
    ax = L0_ctdHelper(meanO2(:,sprIds),"ad",logAxes);
    sgtitle("$O_2$ " + tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/ctd/" + lp + "O2" + tmpT + ".png");

    % summer
    tmpT = "-ad-03";

    chla = load("output\CTD\chla.mat").meanEpN(1:101,sumIds(33:end));
    ax = L0_ctdHelper(chla,"ad",logAxes);
    sgtitle("chl-a " + tmpT);
    exportgraphics(ax,"figures/L0/ctd/"+lp+"chla" + tmpT + ".png");

    % T
    ax = L0_ctdHelper(meanT(:,sumIds),"ad",logAxes);
    sgtitle("T " + tmpT);
    exportgraphics(ax,"figures/L0/ctd/"+lp+"T" + tmpT + ".png");
    
    % Sp
    ax = L0_ctdHelper(meanSp(:,sumIds),"ad",logAxes);
    sgtitle("$S_p$ " + tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/ctd/" + lp + "Sp" + tmpT + ".png");
    
    % O2
    ax = L0_ctdHelper(meanO2(:,sumIds),"ad",logAxes);
    sgtitle("$O_2$ " + tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/ctd/" + lp + "O2" + tmpT + ".png");

    % autumn
    tmpT = "-ad-04";

    chla = load("output\CTD\chla.mat").meanEpN(1:101,sumIds(30:end));
    ax = L0_ctdHelper(chla,"ad",logAxes);
    sgtitle("chl-a " + tmpT);
    exportgraphics(ax,"figures/L0/ctd/"+lp+"chla" + tmpT + ".png");

    % T
    ax = L0_ctdHelper(meanT(:,autIds),"ad",logAxes);
    sgtitle("T " + tmpT);
    exportgraphics(ax,"figures/L0/ctd/"+lp+"T" + tmpT + ".png");
    
    % Sp
    ax = L0_ctdHelper(meanSp(:,autIds),"ad",logAxes);
    sgtitle("$S_p$ " + tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/ctd/" + lp + "Sp" + tmpT + ".png");
    
    % O2
    ax = L0_ctdHelper(meanO2(:,autIds),"ad",logAxes);
    sgtitle("$O_2$ " + tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/ctd/" + lp + "O2" + tmpT + ".png");
end

%% Start Year Sensitivity Analysis.
if startYearAnalysis == true 
    % K-S
    % 2001-2021 (crn 131-)
    tmpT = "";
    chla = load("output\CTD\chla.mat").meanEpN(1:101,131:329);
    ax = L0_ctdHelper(chla);
    sgtitle("chl-a " + tmpT);
    exportgraphics(ax,"figures/L0/ctd/" + lp + "chla" + tmpT + ".png"); clear;
    
    % 2002-2021 (crn 134-)
    tmpT = "02";
    chla = load("output\CTD\chla.mat").meanEpN(1:101,134:329);
    ax = L0_ctdHelper(chla);
    sgtitle("chl-a " + tmpT);
    exportgraphics(ax,"figures/L0/ctd/" + lp + "chla" + tmpT + ".png"); clear;
    
    % 2003-2021 (crn 144-)
    tmpT = "03";
    chla = load("output\CTD\chla.mat").meanEpN(1:101,144:329);
    ax = L0_ctdHelper(chla);
    sgtitle("chl-a " + tmpT);
    exportgraphics(ax,"figures/L0/ctd/" + lp + "chla" + tmpT + ".png"); clear;
    
    % 2004-2021 (crn 155-)
    tmpT = "04";
    chla = load("output\CTD\chla.mat").meanEpN(1:101,155:329);
    ax = L0_ctdHelper(chla);
    sgtitle("chl-a " + tmpT);
    exportgraphics(ax,"figures/L0/ctd/" + lp + "chla" + tmpT + ".png"); clear;
    
    % 2005-2021 (crn 167-)
    tmpT = "05";
    chla = load("output\CTD\chla.mat").meanEpN(1:101,167:329);
    ax = L0_ctdHelper(chla);
    sgtitle("chl-a " + tmpT);
    exportgraphics(ax,"figures/L0/ctd/" + lp + "chla" + tmpT + ".png"); clear;
    
    % 2006-2021 (crn 177-)
    tmpT = "06";
    chla = load("output\CTD\chla.mat").meanEpN(1:101,177:329);
    ax = L0_ctdHelper(chla);
    sgtitle("chl-a " + tmpT);
    exportgraphics(ax,"figures/L0/ctd/" + lp + "chla" + tmpT + ".png"); clear;
    
    % 2007-2021 (crn 189-)
    tmpT = "07";
    chla = load("output\CTD\chla.mat").meanEpN(1:101,189:329);
    ax = L0_ctdHelper(chla);
    sgtitle("chl-a " + tmpT);
    exportgraphics(ax,"figures/L0/ctd/" + lp + "chla" + tmpT + ".png"); clear;
    
    % 2008-2021 (crn 199-)
    tmpT = "08";
    chla = load("output\CTD\chla.mat").meanEpN(1:101,199:329);
    ax = L0_ctdHelper(chla);
    sgtitle("chl-a " + tmpT);
    exportgraphics(ax,"figures/L0/ctd/" + lp + "chla" + tmpT + ".png"); clear;
    
    % 2009-2021 (crn 208-)
    tmpT = "09";
    chla = load("output\CTD\chla.mat").meanEpN(1:101,208:329);
    ax = L0_ctdHelper(chla);
    sgtitle("chl-a " + tmpT);
    exportgraphics(ax,"figures/L0/ctd/" + lp + "chla" + tmpT + ".png"); clear;
    
    % 2010-2021 (crn 219-)
    tmpT = "10";
    chla = load("output\CTD\chla.mat").meanEpN(1:101,219:329);
    ax = L0_ctdHelper(chla);
    sgtitle("chl-a " + tmpT);
    exportgraphics(ax,"figures/L0/ctd/" + lp + "chla" + tmpT + ".png"); clear;
    
    % 2011-2021 (crn 228-)
    tmpT = "11";
    chla = load("output\CTD\chla.mat").meanEpN(1:101,228:329);
    ax = L0_ctdHelper(chla);
    sgtitle("chl-a " + tmpT);
    exportgraphics(ax,"figures/L0/ctd/" + lp + "chla" + tmpT + ".png"); clear;
    
    % 2012-2021 (crn 239-)
    tmpT = "12";
    chla = load("output\CTD\chla.mat").meanEpN(1:101,239:329);
    ax = L0_ctdHelper(chla);
    sgtitle("chl-a " + tmpT);
    exportgraphics(ax,"figures/L0/ctd/" + lp + "chla" + tmpT + ".png"); clear;
    
    % 2013-2021 (crn 249-)
    tmpT = "13";
    chla = load("output\CTD\chla.mat").meanEpN(1:101,249:329);
    ax = L0_ctdHelper(chla);
    sgtitle("chl-a " + tmpT);
    exportgraphics(ax,"figures/L0/ctd/" + lp + "chla" + tmpT + ".png"); clear;
    
    % 2014-2021 (crn 259-)
    tmpT = "14";
    chla = load("output\CTD\chla.mat").meanEpN(1:101,259:329);
    ax = L0_ctdHelper(chla);
    sgtitle("chl-a " + tmpT);
    exportgraphics(ax,"figures/L0/ctd/" + lp + "chla" + tmpT + ".png"); clear;
    
    % 2015-2021 (crn 269-)
    tmpT = "15";
    chla = load("output\CTD\chla.mat").meanEpN(1:101,269:329);
    ax = L0_ctdHelper(chla);
    sgtitle("chl-a " + tmpT);
    exportgraphics(ax,"figures/L0/ctd/" + lp + "chla" + tmpT + ".png"); clear;
    
    % 2016-2021 (crn 280-)
    tmpT = "16";
    chla = load("output\CTD\chla.mat").meanEpN(1:101,280:329);
    ax = L0_ctdHelper(chla);
    sgtitle("chl-a " + tmpT);
    exportgraphics(ax,"figures/L0/ctd/" + lp + "chla" + tmpT + ".png"); clear;
    
    % A-D

    % 2001-2021 (crn 131-)
    tmpT = "ad";
    chla = load("output\CTD\chla.mat").meanEpN(1:101,131:329);
    ax = L0_ctdHelper(chla,"ad");
    sgtitle("chl-a " + tmpT);
    exportgraphics(ax,"figures/L0/ctd/" + lp + "chla" + tmpT + ".png"); clear;
    
    % 2002-2021 (crn 134-)
    tmpT = "02-ad";
    chla = load("output\CTD\chla.mat").meanEpN(1:101,134:329);
    ax = L0_ctdHelper(chla,"ad");
    sgtitle("chl-a " + tmpT);
    exportgraphics(ax,"figures/L0/ctd/" + lp + "chla" + tmpT + ".png"); clear;
    
    % 2003-2021 (crn 144-)
    tmpT = "03-ad";
    chla = load("output\CTD\chla.mat").meanEpN(1:101,144:329);
    ax = L0_ctdHelper(chla,"ad");
    sgtitle("chl-a " + tmpT);
    exportgraphics(ax,"figures/L0/ctd/" + lp + "chla" + tmpT + ".png"); clear;
    
    % 2004-2021 (crn 155-)
    tmpT = "04-ad";
    chla = load("output\CTD\chla.mat").meanEpN(1:101,155:329);
    ax = L0_ctdHelper(chla,"ad");
    sgtitle("chl-a " + tmpT);
    exportgraphics(ax,"figures/L0/ctd/" + lp + "chla" + tmpT + ".png"); clear;
    
    % 2005-2021 (crn 167-)
    tmpT = "05-ad";
    chla = load("output\CTD\chla.mat").meanEpN(1:101,167:329);
    ax = L0_ctdHelper(chla,"ad");
    sgtitle("chl-a " + tmpT);
    exportgraphics(ax,"figures/L0/ctd/" + lp + "chla" + tmpT + ".png"); clear;
    
    % 2006-2021 (crn 177-)
    tmpT = "06-ad";
    chla = load("output\CTD\chla.mat").meanEpN(1:101,177:329);
    ax = L0_ctdHelper(chla,"ad");
    sgtitle("chl-a " + tmpT);
    exportgraphics(ax,"figures/L0/ctd/" + lp + "chla" + tmpT + ".png"); clear;
    
    % 2007-2021 (crn 189-)
    tmpT = "07-ad";
    chla = load("output\CTD\chla.mat").meanEpN(1:101,189:329);
    ax = L0_ctdHelper(chla,"ad");
    sgtitle("chl-a " + tmpT);
    exportgraphics(ax,"figures/L0/ctd/" + lp + "chla" + tmpT + ".png"); clear;
    
    % 2008-2021 (crn 199-)
    tmpT = "08-ad";
    chla = load("output\CTD\chla.mat").meanEpN(1:101,199:329);
    ax = L0_ctdHelper(chla,"ad");
    sgtitle("chl-a " + tmpT);
    exportgraphics(ax,"figures/L0/ctd/" + lp + "chla" + tmpT + ".png"); clear;
    
    % 2009-2021 (crn 208-)
    tmpT = "09-ad";
    chla = load("output\CTD\chla.mat").meanEpN(1:101,208:329);
    ax = L0_ctdHelper(chla,"ad");
    sgtitle("chl-a " + tmpT);
    exportgraphics(ax,"figures/L0/ctd/" + lp + "chla" + tmpT + ".png"); clear;
    
    % 2010-2021 (crn 219-)
    tmpT = "10-ad";
    chla = load("output\CTD\chla.mat").meanEpN(1:101,219:329);
    ax = L0_ctdHelper(chla,"ad");
    sgtitle("chl-a " + tmpT);
    exportgraphics(ax,"figures/L0/ctd/" + lp + "chla" + tmpT + ".png"); clear;
    
    % 2011-2021 (crn 228-)
    tmpT = "11-ad";
    chla = load("output\CTD\chla.mat").meanEpN(1:101,228:329);
    ax = L0_ctdHelper(chla,"ad");
    sgtitle("chl-a " + tmpT);
    exportgraphics(ax,"figures/L0/ctd/" + lp + "chla" + tmpT + ".png"); clear;
    
    % 2012-2021 (crn 239-)
    tmpT = "12-ad";
    chla = load("output\CTD\chla.mat").meanEpN(1:101,239:329);
    ax = L0_ctdHelper(chla,"ad");
    sgtitle("chl-a " + tmpT);
    exportgraphics(ax,"figures/L0/ctd/" + lp + "chla" + tmpT + ".png"); clear;
    
    % 2013-2021 (crn 249-)
    tmpT = "13-ad";
    chla = load("output\CTD\chla.mat").meanEpN(1:101,249:329);
    ax = L0_ctdHelper(chla,"ad");
    sgtitle("chl-a " + tmpT);
    exportgraphics(ax,"figures/L0/ctd/" + lp + "chla" + tmpT + ".png"); clear;
    
    % 2014-2021 (crn 259-)
    tmpT = "14-ad";
    chla = load("output\CTD\chla.mat").meanEpN(1:101,259:329);
    ax = L0_ctdHelper(chla,"ad");
    sgtitle("chl-a " + tmpT);
    exportgraphics(ax,"figures/L0/ctd/" + lp + "chla" + tmpT + ".png"); clear;
    
    % 2015-2021 (crn 269-)
    tmpT = "15-ad";
    chla = load("output\CTD\chla.mat").meanEpN(1:101,269:329);
    ax = L0_ctdHelper(chla,"ad");
    sgtitle("chl-a " + tmpT);
    exportgraphics(ax,"figures/L0/ctd/" + lp + "chla" + tmpT + ".png"); clear;
    
    % 2016-2021 (crn 280-)
    tmpT = "16-ad";
    chla = load("output\CTD\chla.mat").meanEpN(1:101,280:329);
    ax = L0_ctdHelper(chla,"ad");
    sgtitle("chl-a " + tmpT);
    exportgraphics(ax,"figures/L0/ctd/" + lp + "chla" + tmpT + ".png"); clear;

else
    disp("Not running yearly analysis...");
end
%%
% close all;
% 
% tmpAlpha = importdata('data/L0/alphaCaro_50-60.txt');
% tmpBut19 = importdata('data\L0\but19_50-60.txt');
% tmpHex19 = importdata('data\L0\hex19_50-60.txt');
% 
% alphaCar = tmpAlpha.data(:,5);
% but19 = tmpBut19.data(:,5);
% hex19 = tmpHex19.data(:,5);
% 
% %%
% figure; histogram(alphaCar,max(alphaCar)); title('alpha-carotene');
% 
% figure; histogram(but19,max(but19)); title('But-19');
% 
% figure; histogram(hex19,max(hex19)); title('Hex-19');
% 
% %%
% [mleA,ksA] = statsplot2(alphaCar);
% [~,ksB] = statsplot2(but19);
% [~,ksH] = statsplot2(hex19);
% 
% 
% pd = makedist("Gamma","a",mleA(4,1),"b",mleA(4,2));
% 
% [~,pAdA] = adtest(alphaCar,"Distribution",pd);
% [~,pAdB] = adtest(but19,"Distribution",pd);
% [~,pAdH] = adtest(hex19,"Distribution",pd);
% 
% [rA,pA] = bbvuong(alphaCar);
% [rB,pB] = bbvuong(but19);
% [rH,pH] = bbvuong(hex19);
% 
% %%
% dAcar = load("output\L1\acar.mat");
% pOutA = dAcar.pOut;
% cOutA = dAcar.cOut;
% 
% dHex19 = load("output\L1\hex19.mat");
% pOutH = dHex19.pOut;
% cOutH = dHex19.cOut;
% 
% dBut19 = load("output\L1\but19.mat");
% pOutB = dBut19.pOut;
% cOutB = dBut19.cOut;
% %%
% set(groot, 'defaultFigureUnits', 'centimeters', 'defaultFigurePosition', [3 3 36 10]);
% 
% figure;
% subplot(1,3,1)
% histogram(alphaCar); title("alpha-caro");
% subplot(1,3,2)
% histogram(but19); title("but-19");
% subplot(1,3,3)
% histogram(hex19); title("hex-19");
% sgtitle('50-60 dbar bin: L0');
% 
% %%
% figure;
% subplot(1,3,1)
% histogram(cOutA(pOutA==6)); title("alpha-caro");
% subplot(1,3,2)
% histogram(cOutB(pOutB==6)); title("but-19");
% subplot(1,3,3)
% histogram(cOutH(pOutH==6)); title("hex-19");
% sgtitle('50-60 dbar bin: L1');
% 
% %% remove outlier in hex-19 L0
% 
% hex19o = hex19;
% hex19o(40) = [];
% % show
% figure
% plot(hex19);
% hold on
% plot(hex19o);
% hold off
% 
% [~,ksHo] = statsplot2(hex19o);
% [rHo,pHo] = bbvuong(hex19o);