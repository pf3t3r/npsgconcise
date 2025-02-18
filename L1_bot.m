% Script to output L1 bottle results for the statistical analysis.

clear; clc; close all;
addpath("baroneRoutines\");
addpath("func\");

set(groot, "defaultFigureUnits", "centimeters", "defaultFigurePosition", [3 3 15 15]);

% Possible test cases.
principleAnalysisAd = true;     % main analysis: A-D Test
principleAnalysisKs = false;    % main analysis: K-S Test
seasonalAnalysisKs = false;     % seasonality of statistics: K-S
seasonalAnalysisAd = false;     % seasonality of statistics: A-D
testSel = 2; % 2 = norm + logn; 4 = norm + logn + weib + gamm

logAxes = true;                 % output p-values as log values (true)
if logAxes == true
    lp = "log/";
else
    lp = "";
end

%% Extract Maximum Mixed Layer Depth (per cruise) "maxMld"

ctdData = importdata("datafiles\ctd_iso_ALL.mat").ctd;
maxMld = nan(329,1);
for i = 1:329
    if ~isnan([ctdData(i).mld003])
        maxMld(i) = max([ctdData(i).mld003]);
    end
end

clear ctdData i;
save mldVals.mat maxMld;

%% README
% How is this script organised?
% At the highest level we have each variable in "blocks".
% Within these blocks we 
% (1) load the data, (2) extract that data which lies in the mixed layer,
% (3) bin this data, (4) calculate the KS and A-D p-value, skewness, and 
% kurtosis for this binned mixed-layer data, and finally (5) plot (and 
% save?) the results.
% THE TWO SECTIONS BEFORE THE README MUST BE RUN IN ORDER TO LOAD THE MLD,
% THEN EACH INDIVIDUAL BLOCK CAN BE RUN.


%% A-D but with Four Distributions

% A-D
tmpT = "-ad-4";
thresh = 30;

% HPLC chl-a
tmpX = "";
tmp = importdata("data/L1/hplcChla_88-21_150.txt");
ax = L1_helper(tmp,maxMld,thresh,4,"ad");
sgtitle("L1"+tmpX,"Interpreter","latex");
exportgraphics(ax,"figures/L1/bottle/" + lp + "chla" + tmpT + ".png");

%% L1 chla: set up histograms

tmp = importdata("data/L1/hplcChla_88-21_150.txt");
thresh = 30; testSel = 2;
[~,~,~,~,~,~,~,~,~,chla_ML,p_ML] = L1_helper(tmp,maxMld,thresh,testSel,"ad");

% chla = tmp.data(:,5);
% p = tmp.data(:,4);
pB = round(p_ML,-1);   % bin the pressure

%% Plot the L1 chl-a histograms
% for chla measured at pressures in range 5-15 dbar, use pB = 10
% for range 15-25 dbar, use pB = 20, etc.
figure;
histogram(chla_ML(pB==10));

%% histfit variation
figure
histfit(chla_ML(pB==40),10,"lognormal");
%% Principal Analysis: A-D

if principleAnalysisAd == true
    % A-D
    tmpT = "-ad";
    thresh = 30;
   
    % HPLC chl-a
    tmpX = "";
    tmp = importdata("data/L1/hplcChla_88-21_150.txt");
    ax = L1_helper(tmp,maxMld,thresh,testSel,"ad");
    sgtitle("L1"+tmpX,"Interpreter","latex");
    exportgraphics(ax,"figures/L1/bottle/" + lp + "chla" + tmpT + ".png");
    clearvars -except tmpT maxMld lp thresh testSel seasonalAnalysisKs seasonalAnalysisAd;
    
    % mvchla
    tmp = importdata("data/L1/mvchla_94-21_150.txt");
    ax = L1_helper(tmp,maxMld,thresh,testSel,"ad");
    sgtitle("[Monovinyl Chl a] 88-21: L1"+tmpT);
    exportgraphics(ax,"figures/L1/bottle/" + lp + "mvchla" + tmpT + ".png");
    clearvars -except tmpT maxMld lp thresh testSel seasonalAnalysisKs seasonalAnalysisAd;
    
    % dvchla
    tmp = importdata("data/L1/dvchla_94-21_150.txt");
    ax = L1_helper(tmp,maxMld,thresh,testSel,"ad");
    sgtitle("[Divinyl Chl a] 88-21: L1"+tmpT);
    exportgraphics(ax,"figures/L1/bottle/" + lp + "dvchla" + tmpT + ".png");
    clearvars -except tmpT maxMld lp thresh testSel seasonalAnalysisKs seasonalAnalysisAd;
    
    % chl-b
    tmp = importdata("data/L1/chlb_88-21_150.txt");
    ax = L1_helper(tmp,maxMld,thresh,testSel,"ad");
    sgtitle("L1");
    exportgraphics(ax,"figures/L1/bottle/" + lp + "chlb" + tmpT + ".png");
    clearvars -except tmpT maxMld lp thresh testSel seasonalAnalysisKs seasonalAnalysisAd;
    
    % chl-c123
    tmp = importdata("data/L1/chl123_88-21_150.txt");
    ax = L1_helper(tmp,maxMld,thresh,testSel,"ad");
    sgtitle("L1");
    exportgraphics(ax,"figures/L1/bottle/" + lp + "chlc123" + tmpT + ".png");
    clearvars -except tmpT maxMld lp thresh testSel seasonalAnalysisKs seasonalAnalysisAd;
    
    % alpha-caro
    tmp = importdata("data/L1/acar_94-21_150.txt");
    ax = L1_helper(tmp,maxMld,thresh,testSel,"ad");
    sgtitle("alpha-carotene 88-21: L1"+tmpT);
    exportgraphics(ax,"figures/L1/bottle/" + lp + "acar" + tmpT + ".png");
    clearvars -except tmpT maxMld lp thresh testSel seasonalAnalysisKs seasonalAnalysisAd;
    
    % but-19
    tmp = importdata("data/L1/but19_88-21_150.txt");
    ax = L1_helper(tmp,maxMld,thresh,testSel,"ad");
    sgtitle("But-19 88-21: L1"+tmpT);
    exportgraphics(ax,"figures/L1/bottle/" + lp + "but19" + tmpT + ".png");
    clearvars -except tmpT maxMld lp thresh testSel seasonalAnalysisKs seasonalAnalysisAd;
    
    % hex-19
    tmp = importdata("data/L1/hex19_88-21_150.txt");
    ax = L1_helper(tmp,maxMld,thresh,testSel,"ad");
    sgtitle("Hex-19 88-21: L1"+tmpT);
    exportgraphics(ax,"figures/L1/bottle/" + lp + "hex19" + tmpT + ".png");
    clearvars -except tmpT maxMld lp thresh testSel seasonalAnalysisKs seasonalAnalysisAd;
    
    % zeax
    tmp = importdata("data/L1/zeaxan_88-21_150.txt");
    ax = L1_helper(tmp,maxMld,thresh,testSel,"ad");
    sgtitle("Zeaxanthin 88-21: L1"+tmpT);
    exportgraphics(ax,"figures/L1/bottle/" + lp + "zeax" + tmpT + ".png");
    clearvars -except tmpT maxMld lp thresh testSel seasonalAnalysisKs seasonalAnalysisAd;
    
    % pc
    tmpx = "";
    tmp = importdata("data/L1/parc_89-21_150.txt");
    ax = L1_helper(tmp,maxMld,thresh,testSel,"ad");
    sgtitle("L1"+tmpx);
    exportgraphics(ax,"figures/L1/bottle/" + lp + "pc" + tmpT + ".png");
    clearvars -except tmpT maxMld lp thresh testSel seasonalAnalysisKs seasonalAnalysisAd;
    
    % pn
    tmp = importdata("data/L1/parn_89-21_150.txt");
    ax = L1_helper(tmp,maxMld,thresh,testSel,"ad");
    sgtitle("Particulate Nitrogen 88-21: L1"+tmpT);
    exportgraphics(ax,"figures/L1/bottle/" + lp + "pn" + tmpT + ".png");
    clearvars -except tmpT maxMld lp thresh testSel seasonalAnalysisKs seasonalAnalysisAd;
    
    % llp l1
    tmp = importdata("data/L1/llp_88-22_150.txt");
    ax = L1_helper(tmp,maxMld,thresh,testSel,"ad");
    sgtitle("Low-Level Phosphorus 88-21: L1"+tmpT);
    exportgraphics(ax,"figures/L1/bottle/" + lp + "llp" + tmpT + ".png");
    clearvars -except tmpT maxMld lp thresh testSel seasonalAnalysisKs seasonalAnalysisAd;
    
    % lln l1
    tmp = importdata("data/L1/lln_89-22_150.txt");
    ax = L1_helper(tmp,maxMld,thresh,testSel,"ad");
    sgtitle("Low-Level Nitrogen 88-21: L1"+tmpT);
    exportgraphics(ax,"figures/L1/bottle/" + lp + "lln" + tmpT + ".png");
    clearvars -except tmpT maxMld lp thresh testSel seasonalAnalysisKs seasonalAnalysisAd;
    
    % phosphate
    tmp = importdata("data\L1\pho_88-21_150.txt");
    ax = L1_helper(tmp,maxMld,thresh,testSel,"ad");
    sgtitle("Phosphate 88-21: L1"+tmpT);
    exportgraphics(ax,"figures/L1/bottle/" + lp + "phos" + tmpT + ".png");
    clearvars -except tmpT maxMld lp thresh testSel seasonalAnalysisKs seasonalAnalysisAd;
    
    % do
    tmp = importdata("data/L1/oxy_88-21_150.txt");
    ax = L1_helper(tmp,maxMld,thresh,testSel,"ad");
    sgtitle("Dissolved Oxygen 88-21: L1"+tmpT);
    exportgraphics(ax,"figures/L1/bottle/" + lp + "boxy" + tmpT + ".png");
    clearvars -except tmpT maxMld lp thresh testSel seasonalAnalysisKs seasonalAnalysisAd;
    
    % l-12 pp
    tmp = importdata("data/L1/l12_89-22_150.txt");
    ax = L1_helper(tmp,maxMld,thresh,testSel,"ad");
    sgtitle("L-12 PP: L1"+tmpT);
    exportgraphics(ax,"figures/L1/bottle/" + lp + "l12" + tmpT + ".png");
    clearvars -except tmpT maxMld lp thresh testSel seasonalAnalysisKs seasonalAnalysisAd;

end

%% Principal Analysis: K-S

if principleAnalysisKs == true

    % K-S
    thresh = 50;
    tmpT = "";
    
    % chl-a (88-21)
    tmp = importdata("data/L1/hplcChla_88-21_150.txt");
    [ax,p,ks,obs,Sk,Ku,~,~] = L1_helper(tmp,maxMld,thresh,testSel);
    sgtitle("[Chl a] 88-21: L1");
    exportgraphics(ax,"figures/L1/bottle/" + lp + "chla" + tmpT + ".png");
    save("output\L1\chla.mat","p","ks","obs","Sk","Ku");
    clearvars -except tmpT maxMld lp thresh testSel seasonalAnalysisKs seasonalAnalysisAd;
    
    % HPLC Monovinyl chlorophyll a (94-21)
    tmp = importdata("data\L1\mvchla_94-21_150.txt");
    [ax,p,ks,obs,Sk,Ku] = L1_helper(tmp,maxMld,thresh,testSel);
    sgtitle("HPLC Monovinyl chlorophyll a 94-21: L1");
    exportgraphics(ax,"figures/L1/bottle/" + lp + "mvchla" + tmpT + ".png");
    save("output\L1\mvchla.mat","p","ks","obs","Sk","Ku");
    clearvars -except tmpT maxMld lp thresh testSel seasonalAnalysisKs seasonalAnalysisAd;
    
    % HPLC Divinyl chlorophyll a (94-21)
    tmp = importdata("data\L1\dvchla_94-21_150.txt");
    [ax,p,ks,obs,Sk,Ku] = L1_helper(tmp,maxMld,thresh,testSel);
    sgtitle("HPLC Divinyl chlorophyll a 94-21: L1");
    exportgraphics(ax,"figures/L1/bottle/" + lp + "dvchla" + tmpT + ".png");
    save("output\L1\dvchla.mat","p","ks","obs","Sk","Ku");
    clearvars -except tmpT maxMld lp thresh testSel seasonalAnalysisKs seasonalAnalysisAd;
    
    % HPLC chlorophyll b (88-21)
    tmp = importdata("data\L1\chlb_88-21_150.txt");
    [ax,p,ks,obs,Sk,Ku] = L1_helper(tmp,maxMld,thresh,testSel);
    sgtitle("HPLC chlorophyll b 88-21: L1");
    exportgraphics(ax,"figures/L1/bottle/" + lp + "chlb" + tmpT + ".png");
    save("output\L1\chlb.mat","p","ks","obs","Sk","Ku");
    clearvars -except tmpT maxMld lp thresh testSel seasonalAnalysisKs seasonalAnalysisAd;
    
    % HPLC Chlorophyll C1 + C2 + C3 (88-21)
    tmp = importdata("data\L1\chl123_88-21_150.txt");
    [ax,p,ks,obs,Sk,Ku] = L1_helper(tmp,maxMld,thresh,testSel);
    sgtitle("HPLC Chl c1 + c2 + c3: L1");
    exportgraphics(ax,"figures/L1/bottle/" + lp + "chlc123" + tmpT + ".png");
    save("output\L1\chl123.mat","p","ks","obs","Sk","Ku");
    clearvars -except tmpT maxMld lp thresh testSel seasonalAnalysisKs seasonalAnalysisAd;
    
    % HPLC alpha-Carotene (94-21)
    tmp = importdata("data\L1\acar_94-21_150.txt");
    [ax,p,ks,obs,Sk,Ku] = L1_helper(tmp,maxMld,thresh,testSel);
    sgtitle("HPLC alpha-Carotene 94-21: L1");
    exportgraphics(ax,"figures/L1/bottle/" + lp + "acar" + tmpT + ".png");
    save("output\L1\acar.mat","p","ks","obs","Sk","Ku");
    clearvars -except tmpT maxMld lp thresh testSel seasonalAnalysisKs seasonalAnalysisAd;
    
    % HPLC 19' Butanoyloxyfucoxanthin (88-21)
    tmp = importdata("data\L1\but19_88-21_150.txt");
    [ax,p,ks,obs,Sk,Ku,~,~,~] = L1_helper(tmp,maxMld,thresh,testSel);
    sgtitle("HPLC 19 Butanoyloxyfucoxanthin 88-21: L1");
    exportgraphics(ax,"figures/L1/bottle/" + lp + "but19" + tmpT + ".png");
    save("output\L1\but19.mat","p","ks","obs","Sk","Ku");
    clearvars -except tmpT maxMld lp thresh testSel seasonalAnalysisKs seasonalAnalysisAd;
    
    % HPLC 19" Hexanoyloxyfucoxanthin (88-21)
    tmp = importdata("data\L1\hex19_88-21_150.txt");
    [ax,p,ks,obs,Sk,Ku,~,~,~] = L1_helper(tmp,maxMld,thresh,testSel);
    sgtitle("HPLC 19' Hexanoyloxyfucoxanthin 88-21: L1");
    exportgraphics(ax,"figures/L1/bottle/" + lp + "hex19" + tmpT + ".png");
    save("output\L1\hex19.mat","p","ks","obs","Sk","Ku");
    clearvars -except tmpT maxMld lp thresh testSel seasonalAnalysisKs seasonalAnalysisAd;
    
    % HPLC Zeaxanthin (88-21)
    tmp = importdata("data\L1\zeaxan_88-21_150.txt");
    [ax,p,ks,obs,Sk,Ku] = L1_helper(tmp,maxMld,thresh,testSel);
    sgtitle("HPLC Zeaxanthin 88-21: L1");
    exportgraphics(ax,"figures/L1/bottle/" + lp + "zeax" + tmpT + ".png");
    save("output\L1\zeaxan.mat","p","ks","obs","Sk","Ku");
    clearvars -except tmpT maxMld lp thresh testSel seasonalAnalysisKs seasonalAnalysisAd;
    
    % Particulate Carbon (89-21)
    tmp = importdata("data\L1\parc_89-21_150.txt");
    [ax,p,ks,obs,Sk,Ku] = L1_helper(tmp,maxMld,thresh,testSel);
    sgtitle("Particulate Carbon 89-21: L1");
    exportgraphics(ax,"figures/L1/bottle/" + lp + "pc" + tmpT + ".png");
    save("output\L1\pc.mat","p","ks","obs","Sk","Ku");
    clearvars -except tmpT maxMld lp thresh testSel seasonalAnalysisKs seasonalAnalysisAd;
    
    % Particulate Nitrogen (89-21)
    tmp = importdata("data\L1\parn_89-21_150.txt");
    [ax,p,ks,obs,Sk,Ku] = L1_helper(tmp,maxMld,thresh,testSel);
    sgtitle("Particulate Nitrogen 89-21: L1");
    exportgraphics(ax,"figures/L1/bottle/" + lp + "pn" + tmpT + ".png");
    save("output\L1\pn.mat","p","ks","obs","Sk","Ku");
    clearvars -except tmpT maxMld lp thresh testSel seasonalAnalysisKs seasonalAnalysisAd;
    
    % Low-Level Phosphorus (88-22)
    tmp = importdata("data\L1\llp_88-22_150.txt");
    [ax,p,ks,obs,Sk,Ku] = L1_helper(tmp,maxMld,thresh,testSel);
    sgtitle("Low-Level Phosphorus 88-22: L1");
    exportgraphics(ax,"figures/L1/bottle/" + lp + "llp" + tmpT + ".png");
    save("output\L1\llp.mat","p","ks","obs","Sk","Ku");
    clearvars -except tmpT maxMld lp thresh testSel seasonalAnalysisKs seasonalAnalysisAd;
    
    % Low-Level Nitrogen (89-22)
    tmp = importdata("data\L1\lln_89-22_150.txt");
    [ax,p,ks,obs,Sk,Ku] = L1_helper(tmp,maxMld,thresh,testSel);
    sgtitle("Low-Level Nitrogen 89-22: L1");
    exportgraphics(ax,"figures/L1/bottle/" + lp + "lln" + tmpT + ".png");
    save("output\L1\lln.mat","p","ks","obs","Sk","Ku");
    clearvars -except tmpT maxMld lp thresh testSel seasonalAnalysisKs seasonalAnalysisAd;
    
    % Bottle Dissolved Oxygen (88-21)
    tmp = importdata("data/L1/oxy_88-21_150.txt");
    [ax,p,ks,obs,Sk,Ku] = L1_helper(tmp,maxMld,thresh,testSel);
    sgtitle("Dissolved Oxygen 88-21: L1");
    exportgraphics(ax,"figures/L1/bottle/" + lp + "boxy" + tmpT + ".png");
    save("output\L1\boxy.mat","p","ks","obs","Sk","Ku");
    clearvars -except tmpT maxMld lp thresh testSel seasonalAnalysisKs seasonalAnalysisAd;
    
    % PProd Light-12 (89-22) (4 significant digits => very good!)
    tmp = importdata("data\L1\l12_89-22_150.txt");
    [ax,p,ks,obs,Sk,Ku] = L1_helper(tmp,maxMld,thresh,testSel);
    sgtitle("PProd Light-12 89-22: L1");
    exportgraphics(ax,"figures/L1/bottle/" + lp + "l12" + tmpT + ".png");
    save("output\L1\l12.mat","p","ks","obs","Sk","Ku");
    clearvars -except tmpT maxMld lp thresh testSel seasonalAnalysisKs seasonalAnalysisAd; 
end

%% Seasonal Analysis: K-S
if seasonalAnalysisKs == true
    
    thresh = 50;
    % WINTER
    tmpT = "-01";

    % chla
    tmp = importdata("data/L1/hplcChla_88-21_150.txt");
    ax = L1_helper(tmp,maxMld,thresh,testSel,"ks",true,1);
    sgtitle("[Chl a] 88-21: L1" + tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L1/bottle/" + lp + "chla" + tmpT + ".png");
    clear tmp ax;

    % mvchla (94-21)
    tmp = importdata("data\L1\mvchla_94-21_150.txt");
    ax = L1_helper(tmp,maxMld,thresh,testSel,"ks",true,1);
    sgtitle("HPLC Monovinyl chlorophyll a 94-21: L1"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L1/bottle/" + lp + "mvchla" + tmpT + ".png");
    clear tmp ax;
    
    % dvchla (94-21)
    tmp = importdata("data\L1\dvchla_94-21_150.txt");
    ax = L1_helper(tmp,maxMld,thresh,testSel,"ks",true,1);
    sgtitle("HPLC Divinyl chlorophyll a 94-21: L1"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L1/bottle/" + lp + "dvchla" + tmpT + ".png");
    clear tmp ax;
    
    % chlb (88-21)
    tmp = importdata("data\L1\chlb_88-21_150.txt");
    ax = L1_helper(tmp,maxMld,thresh,testSel,"ks",true,1);
    sgtitle("HPLC chlorophyll b 88-21: L1"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L1/bottle/" + lp + "chlb" + tmpT + ".png");
    clear tmp ax;
    
    % HPLC Chlorophyll C1 + C2 + C3 (88-21)
    tmp = importdata("data\L1\chl123_88-21_150.txt");
    ax = L1_helper(tmp,maxMld,thresh,testSel,"ks",true,1);
    sgtitle("HPLC Chl c1 + c2 + c3: L1"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L1/bottle/" + lp + "chlc123" + tmpT + ".png");
    clear tmp ax;
    
    % HPLC alpha-Carotene (94-21)
    tmp = importdata("data\L1\acar_94-21_150.txt");
    ax = L1_helper(tmp,maxMld,thresh,testSel,"ks",true,1);
    sgtitle("HPLC alpha-Carotene 94-21: L1"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L1/bottle/" + lp + "acar" + tmpT + ".png");
    clear tmp ax;
    
    % HPLC 19' Butanoyloxyfucoxanthin (88-21)
    tmp = importdata("data\L1\but19_88-21_150.txt");
    ax = L1_helper(tmp,maxMld,thresh,testSel,"ks",true,1);
    sgtitle("HPLC 19 Butanoyloxyfucoxanthin 88-21: L1"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L1/bottle/" + lp + "but19" + tmpT + ".png");
    clear tmp ax;
    
    % HPLC 19" Hexanoyloxyfucoxanthin (88-21)
    tmp = importdata("data\L1\hex19_88-21_150.txt");
    ax = L1_helper(tmp,maxMld,thresh,testSel,"ks",true,1);
    sgtitle("HPLC 19' Hexanoyloxyfucoxanthin 88-21: L1"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L1/bottle/" + lp + "hex19" + tmpT + ".png");
    clear tmp ax;
    
    % HPLC Zeaxanthin (88-21)
    tmp = importdata("data\L1\zeaxan_88-21_150.txt");
    ax = L1_helper(tmp,maxMld,thresh,testSel,"ks",true,1);
    sgtitle("HPLC Zeaxanthin 88-21: L1"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L1/bottle/" + lp + "zeax" + tmpT + ".png");
    clear tmp ax;
    
    % Particulate Carbon (89-21)
    tmp = importdata("data\L1\parc_89-21_150.txt");
    ax = L1_helper(tmp,maxMld,thresh,testSel,"ks",true,1);
    sgtitle("Particulate Carbon 89-21: L1"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L1/bottle/" + lp + "pc" + tmpT + ".png");
    clear tmp ax;
    
    % Particulate Nitrogen (89-21)
    tmp = importdata("data\L1\parn_89-21_150.txt");
    ax = L1_helper(tmp,maxMld,thresh,testSel,"ks",true,1);
    sgtitle("Particulate Nitrogen 89-21: L1"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L1/bottle/" + lp + "pn" + tmpT + ".png");
    clear tmp ax;
    
    % Low-Level Phosphorus (88-22)
    tmp = importdata("data\L1\llp_88-22_150.txt");
    ax = L1_helper(tmp,maxMld,thresh,testSel,"ks",true,1);
    sgtitle("Low-Level Phosphorus 88-22: L1"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L1/bottle/" + lp + "llp" + tmpT + ".png");
    clear tmp ax;
    
    % Low-Level Nitrogen (89-22)
    tmp = importdata("data\L1\lln_89-22_150.txt");
    ax = L1_helper(tmp,maxMld,thresh,testSel,"ks",true,1);
    sgtitle("Low-Level Nitrogen 89-22: L1"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L1/bottle/" + lp + "lln" + tmpT + ".png");
    clear tmp ax;
    
    % Bottle Dissolved Oxygen (88-21)
    tmp = importdata("data/L1/oxy_88-21_150.txt");
    ax = L1_helper(tmp,maxMld,thresh,testSel,"ks",true,1);
    sgtitle("Dissolved Oxygen 88-21: L1"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L1/bottle/" + lp + "boxy" + tmpT + ".png");
    clear tmp ax;
    
    % PProd Light-12 (89-22) (4 significant digits => very good!)
    tmp = importdata("data\L1\l12_89-22_150.txt");
    ax = L1_helper(tmp,maxMld,thresh,testSel,"ks",true,1);
    sgtitle("PProd Light-12 89-22: L1"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L1/bottle/" + lp + "l12" + tmpT + ".png");
    clear tmp ax;

    %%% SPRING
    tmpT = "-02";

    % chla
    tmp = importdata("data/L1/hplcChla_88-21_150.txt");
    ax = L1_helper(tmp,maxMld,thresh,testSel,"ks",true,2);
    sgtitle("[Chl a] 88-21: L1" + tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L1/bottle/" + lp + "chla" + tmpT + ".png");
    clear tmp ax;

    % mvchla (94-21)
    tmp = importdata("data\L1\mvchla_94-21_150.txt");
    ax = L1_helper(tmp,maxMld,thresh,testSel,"ks",true,2);
    sgtitle("HPLC Monovinyl chlorophyll a 94-21: L1"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L1/bottle/" + lp + "mvchla" + tmpT + ".png");
    clear tmp ax;
    
    % dvchla (94-21)
    tmp = importdata("data\L1\dvchla_94-21_150.txt");
    ax = L1_helper(tmp,maxMld,thresh,testSel,"ks",true,2);
    sgtitle("HPLC Divinyl chlorophyll a 94-21: L1"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L1/bottle/" + lp + "dvchla" + tmpT + ".png");
    clear tmp ax;
    
    % chlb (88-21)
    tmp = importdata("data\L1\chlb_88-21_150.txt");
    ax = L1_helper(tmp,maxMld,thresh,testSel,"ks",true,2);
    sgtitle("HPLC chlorophyll b 88-21: L1"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L1/bottle/" + lp + "chlb" + tmpT + ".png");
    clear tmp ax;
    
    % HPLC Chlorophyll C1 + C2 + C3 (88-21)
    tmp = importdata("data\L1\chl123_88-21_150.txt");
    ax = L1_helper(tmp,maxMld,thresh,testSel,"ks",true,2);
    sgtitle("HPLC Chl c1 + c2 + c3: L1"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L1/bottle/" + lp + "chlc123" + tmpT + ".png");
    clear tmp ax;
    
    % HPLC alpha-Carotene (94-21)
    tmp = importdata("data\L1\acar_94-21_150.txt");
    ax = L1_helper(tmp,maxMld,thresh,testSel,"ks",true,2);
    sgtitle("HPLC alpha-Carotene 94-21: L1"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L1/bottle/" + lp + "acar" + tmpT + ".png");
    clear tmp ax;
    
    % HPLC 19' Butanoyloxyfucoxanthin (88-21)
    tmp = importdata("data\L1\but19_88-21_150.txt");
    ax = L1_helper(tmp,maxMld,thresh,testSel,"ks",true,2);
    sgtitle("HPLC 19 Butanoyloxyfucoxanthin 88-21: L1"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L1/bottle/" + lp + "but19" + tmpT + ".png");
    clear tmp ax;
    
    % HPLC 19" Hexanoyloxyfucoxanthin (88-21)
    tmp = importdata("data\L1\hex19_88-21_150.txt");
    ax = L1_helper(tmp,maxMld,thresh,testSel,"ks",true,2);
    sgtitle("HPLC 19' Hexanoyloxyfucoxanthin 88-21: L1"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L1/bottle/" + lp + "hex19" + tmpT + ".png");
    clear tmp ax;
    
    % HPLC Zeaxanthin (88-21)
    tmp = importdata("data\L1\zeaxan_88-21_150.txt");
    ax = L1_helper(tmp,maxMld,thresh,testSel,"ks",true,2);
    sgtitle("HPLC Zeaxanthin 88-21: L1"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L1/bottle/" + lp + "zeax" + tmpT + ".png");
    clear tmp ax;
    
    % Particulate Carbon (89-21)
    tmp = importdata("data\L1\parc_89-21_150.txt");
    ax = L1_helper(tmp,maxMld,thresh,testSel,"ks",true,2);
    sgtitle("Particulate Carbon 89-21: L1"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L1/bottle/" + lp + "pc" + tmpT + ".png");
    clear tmp ax;
    
    % Particulate Nitrogen (89-21)
    tmp = importdata("data\L1\parn_89-21_150.txt");
    ax = L1_helper(tmp,maxMld,thresh,testSel,"ks",true,2);
    sgtitle("Particulate Nitrogen 89-21: L1"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L1/bottle/" + lp + "pn" + tmpT + ".png");
    clear tmp ax;
    
    % Low-Level Phosphorus (88-22)
    tmp = importdata("data\L1\llp_88-22_150.txt");
    ax = L1_helper(tmp,maxMld,thresh,testSel,"ks",true,2);
    sgtitle("Low-Level Phosphorus 88-22: L1"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L1/bottle/" + lp + "llp" + tmpT + ".png");
    clear tmp ax;
    
    % Low-Level Nitrogen (89-22)
    tmp = importdata("data\L1\lln_89-22_150.txt");
    ax = L1_helper(tmp,maxMld,thresh,testSel,"ks",true,2);
    sgtitle("Low-Level Nitrogen 89-22: L1"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L1/bottle/" + lp + "lln" + tmpT + ".png");
    clear tmp ax;
    
    % Bottle Dissolved Oxygen (88-21)
    tmp = importdata("data/L1/oxy_88-21_150.txt");
    ax = L1_helper(tmp,maxMld,thresh,testSel,"ks",true,2);
    sgtitle("Dissolved Oxygen 88-21: L1"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L1/bottle/" + lp + "boxy" + tmpT + ".png");
    clear tmp ax;
    
    % PProd Light-12 (89-22) (4 significant digits => very good!)
    tmp = importdata("data\L1\l12_89-22_150.txt");
    ax = L1_helper(tmp,maxMld,thresh,testSel,"ks",true,2);
    sgtitle("PProd Light-12 89-22: L1"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L1/bottle/" + lp + "l12" + tmpT + ".png");
    clear tmp ax;

    % summer
    tmpT = "-03";

    % chla
    tmp = importdata("data/L1/hplcChla_88-21_150.txt");
    ax = L1_helper(tmp,maxMld,thresh,testSel,"ks",true,3);
    sgtitle("[Chl a] 88-21: L1" + tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L1/bottle/" + lp + "chla" + tmpT + ".png");
    clear tmp ax;

    % mvchla (94-21)
    tmp = importdata("data\L1\mvchla_94-21_150.txt");
    ax = L1_helper(tmp,maxMld,thresh,testSel,"ks",true,3);
    sgtitle("HPLC Monovinyl chlorophyll a 94-21: L1"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L1/bottle/" + lp + "mvchla" + tmpT + ".png");
    clear tmp ax;
    
    % dvchla (94-21)
    tmp = importdata("data\L1\dvchla_94-21_150.txt");
    ax = L1_helper(tmp,maxMld,thresh,testSel,"ks",true,3);
    sgtitle("HPLC Divinyl chlorophyll a 94-21: L1"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L1/bottle/" + lp + "dvchla" + tmpT + ".png");
    clear tmp ax;
    
    % chlb (88-21)
    tmp = importdata("data\L1\chlb_88-21_150.txt");
    ax = L1_helper(tmp,maxMld,thresh,testSel,"ks",true,3);
    sgtitle("HPLC chlorophyll b 88-21: L1"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L1/bottle/" + lp + "chlb" + tmpT + ".png");
    clear tmp ax;
    
    % HPLC Chlorophyll C1 + C2 + C3 (88-21)
    tmp = importdata("data\L1\chl123_88-21_150.txt");
    ax = L1_helper(tmp,maxMld,thresh,testSel,"ks",true,3);
    sgtitle("HPLC Chl c1 + c2 + c3: L1"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L1/bottle/" + lp + "chlc123" + tmpT + ".png");
    clear tmp ax;
    
    % HPLC alpha-Carotene (94-21)
    tmp = importdata("data\L1\acar_94-21_150.txt");
    ax = L1_helper(tmp,maxMld,thresh,testSel,"ks",true,3);
    sgtitle("HPLC alpha-Carotene 94-21: L1"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L1/bottle/" + lp + "acar" + tmpT + ".png");
    clear tmp ax;
    
    % HPLC 19' Butanoyloxyfucoxanthin (88-21)
    tmp = importdata("data\L1\but19_88-21_150.txt");
    ax = L1_helper(tmp,maxMld,thresh,testSel,"ks",true,3);
    sgtitle("HPLC 19 Butanoyloxyfucoxanthin 88-21: L1"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L1/bottle/" + lp + "but19" + tmpT + ".png");
    clear tmp ax;
    
    % HPLC 19" Hexanoyloxyfucoxanthin (88-21)
    tmp = importdata("data\L1\hex19_88-21_150.txt");
    ax = L1_helper(tmp,maxMld,thresh,testSel,"ks",true,3);
    sgtitle("HPLC 19' Hexanoyloxyfucoxanthin 88-21: L1"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L1/bottle/" + lp + "hex19" + tmpT + ".png");
    clear tmp ax;
    
    % HPLC Zeaxanthin (88-21)
    tmp = importdata("data\L1\zeaxan_88-21_150.txt");
    ax = L1_helper(tmp,maxMld,thresh,testSel,"ks",true,3);
    sgtitle("HPLC Zeaxanthin 88-21: L1"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L1/bottle/" + lp + "zeax" + tmpT + ".png");
    clear tmp ax;
    
    % Particulate Carbon (89-21)
    tmp = importdata("data\L1\parc_89-21_150.txt");
    ax = L1_helper(tmp,maxMld,thresh,testSel,"ks",true,3);
    sgtitle("Particulate Carbon 89-21: L1"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L1/bottle/" + lp + "pc" + tmpT + ".png");
    clear tmp ax;
    
    % Particulate Nitrogen (89-21)
    tmp = importdata("data\L1\parn_89-21_150.txt");
    ax = L1_helper(tmp,maxMld,thresh,testSel,"ks",true,3);
    sgtitle("Particulate Nitrogen 89-21: L1"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L1/bottle/" + lp + "pn" + tmpT + ".png");
    clear tmp ax;
    
    % Low-Level Phosphorus (88-22)
    tmp = importdata("data\L1\llp_88-22_150.txt");
    ax = L1_helper(tmp,maxMld,thresh,testSel,"ks",true,3);
    sgtitle("Low-Level Phosphorus 88-22: L1"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L1/bottle/" + lp + "llp" + tmpT + ".png");
    clear tmp ax;
    
    % Low-Level Nitrogen (89-22)
    tmp = importdata("data\L1\lln_89-22_150.txt");
    ax = L1_helper(tmp,maxMld,thresh,testSel,"ks",true,3);
    sgtitle("Low-Level Nitrogen 89-22: L1"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L1/bottle/" + lp + "lln" + tmpT + ".png");
    clear tmp ax;
    
    % Bottle Dissolved Oxygen (88-21)
    tmp = importdata("data/L1/oxy_88-21_150.txt");
    ax = L1_helper(tmp,maxMld,thresh,testSel,"ks",true,3);
    sgtitle("Dissolved Oxygen 88-21: L1"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L1/bottle/" + lp + "boxy" + tmpT + ".png");
    clear tmp ax;
    
    % PProd Light-12 (89-22) (4 significant digits => very good!)
    tmp = importdata("data\L1\l12_89-22_150.txt");
    ax = L1_helper(tmp,maxMld,thresh,testSel,"ks",true,3);
    sgtitle("PProd Light-12 89-22: L1"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L1/bottle/" + lp + "l12" + tmpT + ".png");
    clear tmp ax;

    % autumn
    tmpT = "-04";

    % chla
    tmp = importdata("data/L1/hplcChla_88-21_150.txt");
    ax = L1_helper(tmp,maxMld,thresh,testSel,"ks",true,4);
    sgtitle("[Chl a] 88-21: L1" + tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L1/bottle/" + lp + "chla" + tmpT + ".png");
    clear tmp ax;

    % mvchla (94-21)
    tmp = importdata("data\L1\mvchla_94-21_150.txt");
    ax = L1_helper(tmp,maxMld,thresh,testSel,"ks",true,4);
    sgtitle("HPLC Monovinyl chlorophyll a 94-21: L1"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L1/bottle/" + lp + "mvchla" + tmpT + ".png");
    clear tmp ax;
    
    % dvchla (94-21)
    tmp = importdata("data\L1\dvchla_94-21_150.txt");
    ax = L1_helper(tmp,maxMld,thresh,testSel,"ks",true,4);
    sgtitle("HPLC Divinyl chlorophyll a 94-21: L1"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L1/bottle/" + lp + "dvchla" + tmpT + ".png");
    clear tmp ax;
    
    % chlb (88-21)
    tmp = importdata("data\L1\chlb_88-21_150.txt");
    ax = L1_helper(tmp,maxMld,thresh,testSel,"ks",true,4);
    sgtitle("HPLC chlorophyll b 88-21: L1"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L1/bottle/" + lp + "chlb" + tmpT + ".png");
    clear tmp ax;
    
    % HPLC Chlorophyll C1 + C2 + C3 (88-21)
    tmp = importdata("data\L1\chl123_88-21_150.txt");
    ax = L1_helper(tmp,maxMld,thresh,testSel,"ks",true,4);
    sgtitle("HPLC Chl c1 + c2 + c3: L1"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L1/bottle/" + lp + "chlc123" + tmpT + ".png");
    clear tmp ax;
    
    % HPLC alpha-Carotene (94-21)
    tmp = importdata("data\L1\acar_94-21_150.txt");
    ax = L1_helper(tmp,maxMld,thresh,testSel,"ks",true,4);
    sgtitle("HPLC alpha-Carotene 94-21: L1"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L1/bottle/" + lp + "acar" + tmpT + ".png");
    clear tmp ax;
    
    % HPLC 19' Butanoyloxyfucoxanthin (88-21)
    tmp = importdata("data\L1\but19_88-21_150.txt");
    ax = L1_helper(tmp,maxMld,thresh,testSel,"ks",true,4);
    sgtitle("HPLC 19 Butanoyloxyfucoxanthin 88-21: L1"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L1/bottle/" + lp + "but19" + tmpT + ".png");
    clear tmp ax;
    
    % HPLC 19" Hexanoyloxyfucoxanthin (88-21)
    tmp = importdata("data\L1\hex19_88-21_150.txt");
    ax = L1_helper(tmp,maxMld,thresh,testSel,"ks",true,4);
    sgtitle("HPLC 19' Hexanoyloxyfucoxanthin 88-21: L1"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L1/bottle/" + lp + "hex19" + tmpT + ".png");
    clear tmp ax;
    
    % HPLC Zeaxanthin (88-21)
    tmp = importdata("data\L1\zeaxan_88-21_150.txt");
    ax = L1_helper(tmp,maxMld,thresh,testSel,"ks",true,4);
    sgtitle("HPLC Zeaxanthin 88-21: L1"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L1/bottle/" + lp + "zeax" + tmpT + ".png");
    clear tmp ax;
    
    % Particulate Carbon (89-21)
    tmp = importdata("data\L1\parc_89-21_150.txt");
    ax = L1_helper(tmp,maxMld,thresh,testSel,"ks",true,4);
    sgtitle("Particulate Carbon 89-21: L1"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L1/bottle/" + lp + "pc" + tmpT + ".png");
    clear tmp ax;
    
    % Particulate Nitrogen (89-21)
    tmp = importdata("data\L1\parn_89-21_150.txt");
    ax = L1_helper(tmp,maxMld,thresh,testSel,"ks",true,4);
    sgtitle("Particulate Nitrogen 89-21: L1"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L1/bottle/" + lp + "pn" + tmpT + ".png");
    clear tmp ax;
    
    % Low-Level Phosphorus (88-22)
    tmp = importdata("data\L1\llp_88-22_150.txt");
    ax = L1_helper(tmp,maxMld,thresh,testSel,"ks",true,4);
    sgtitle("Low-Level Phosphorus 88-22: L1"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L1/bottle/" + lp + "llp" + tmpT + ".png");
    clear tmp ax;
    
    % Low-Level Nitrogen (89-22)
    tmp = importdata("data\L1\lln_89-22_150.txt");
    ax = L1_helper(tmp,maxMld,thresh,testSel,"ks",true,4);
    sgtitle("Low-Level Nitrogen 89-22: L1"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L1/bottle/" + lp + "lln" + tmpT + ".png");
    clear tmp ax;
    
    % Bottle Dissolved Oxygen (88-21)
    tmp = importdata("data/L1/oxy_88-21_150.txt");
    ax = L1_helper(tmp,maxMld,thresh,testSel,"ks",true,4);
    sgtitle("Dissolved Oxygen 88-21: L1"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L1/bottle/" + lp + "boxy" + tmpT + ".png");
    clear tmp ax;
    
    % PProd Light-12 (89-22) (4 significant digits => very good!)
    tmp = importdata("data\L1\l12_89-22_150.txt");
    ax = L1_helper(tmp,maxMld,thresh,testSel,"ks",true,4);
    sgtitle("PProd Light-12 89-22: L1"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L1/bottle/" + lp + "l12" + tmpT + ".png");
    clear tmp ax;

end

%% Seasonal Analysis: A-D

if seasonalAnalysisAd == true

    thresh = 30;
    set(groot, "defaultFigureUnits", "centimeters", "defaultFigurePosition", [3 3 10 15]);

    %%% WINTER
    tmpT = "-ad-01";

    % chla
    tmpx = "Winter";
    tmp = importdata("data/L1/hplcChla_88-21_150.txt");
    ax = L1_helper(tmp,maxMld,thresh,testSel,"ad",true,1);
    sgtitle("L1 " + tmpx,"Interpreter","latex");
    exportgraphics(ax,"figures/L1/bottle/" + lp + "chla" + tmpT + ".png");
    clear tmp ax;

    % mvchla (94-21)
    tmp = importdata("data\L1\mvchla_94-21_150.txt");
    ax = L1_helper(tmp,maxMld,thresh,testSel,"ad",true,1);
    sgtitle("HPLC Monovinyl chlorophyll a 94-21: L1"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L1/bottle/" + lp + "mvchla" + tmpT + ".png");
    clear tmp ax;
    
    % dvchla (94-21)
    tmp = importdata("data\L1\dvchla_94-21_150.txt");
    ax = L1_helper(tmp,maxMld,thresh,testSel,"ad",true,1);
    sgtitle("HPLC Divinyl chlorophyll a 94-21: L1"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L1/bottle/" + lp + "dvchla" + tmpT + ".png");
    clear tmp ax;
    
    % chlb (88-21)
    tmp = importdata("data\L1\chlb_88-21_150.txt");
    ax = L1_helper(tmp,maxMld,thresh,testSel,"ad",true,1);
    sgtitle("HPLC chlorophyll b 88-21: L1"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L1/bottle/" + lp + "chlb" + tmpT + ".png");
    clear tmp ax;
    
    % HPLC Chlorophyll C1 + C2 + C3 (88-21)
    tmp = importdata("data\L1\chl123_88-21_150.txt");
    ax = L1_helper(tmp,maxMld,thresh,testSel,"ad",true,1);
    sgtitle("HPLC Chl c1 + c2 + c3: L1"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L1/bottle/" + lp + "chlc123" + tmpT + ".png");
    clear tmp ax;
    
    % HPLC alpha-Carotene (94-21)
    tmp = importdata("data\L1\acar_94-21_150.txt");
    ax = L1_helper(tmp,maxMld,thresh,testSel,"ad",true,1);
    sgtitle("HPLC alpha-Carotene 94-21: L1"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L1/bottle/" + lp + "acar" + tmpT + ".png");
    clear tmp ax;
    
    % HPLC 19' Butanoyloxyfucoxanthin (88-21)
    tmp = importdata("data\L1\but19_88-21_150.txt");
    ax = L1_helper(tmp,maxMld,thresh,testSel,"ad",true,1);
    sgtitle("HPLC 19 Butanoyloxyfucoxanthin 88-21: L1"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L1/bottle/" + lp + "but19" + tmpT + ".png");
    clear tmp ax;
    
    % HPLC 19" Hexanoyloxyfucoxanthin (88-21)
    tmp = importdata("data\L1\hex19_88-21_150.txt");
    ax = L1_helper(tmp,maxMld,thresh,testSel,"ad",true,1);
    sgtitle("HPLC 19' Hexanoyloxyfucoxanthin 88-21: L1"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L1/bottle/" + lp + "hex19" + tmpT + ".png");
    clear tmp ax;
    
    % HPLC Zeaxanthin (88-21)
    tmp = importdata("data\L1\zeaxan_88-21_150.txt");
    ax = L1_helper(tmp,maxMld,thresh,testSel,"ad",true,1);
    sgtitle("HPLC Zeaxanthin 88-21: L1"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L1/bottle/" + lp + "zeax" + tmpT + ".png");
    clear tmp ax;
    
    % Particulate Carbon (89-21)
    tmp = importdata("data\L1\parc_89-21_150.txt");
    ax = L1_helper(tmp,maxMld,thresh,testSel,"ad",true,1);
    sgtitle("L1 Winter","Interpreter","latex");
    exportgraphics(ax,"figures/L1/bottle/" + lp + "pc" + tmpT + ".png");
    clear tmp ax;
    
    % Particulate Nitrogen (89-21)
    tmp = importdata("data\L1\parn_89-21_150.txt");
    ax = L1_helper(tmp,maxMld,thresh,testSel,"ad",true,1);
    sgtitle("Particulate Nitrogen 89-21: L1"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L1/bottle/" + lp + "pn" + tmpT + ".png");
    clear tmp ax;
    
    % Low-Level Phosphorus (88-22)
    tmp = importdata("data\L1\llp_88-22_150.txt");
    ax = L1_helper(tmp,maxMld,thresh,testSel,"ad",true,1);
    sgtitle("Low-Level Phosphorus 88-22: L1"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L1/bottle/" + lp + "llp" + tmpT + ".png");
    clear tmp ax;
    
    % Low-Level Nitrogen (89-22)
    tmp = importdata("data\L1\lln_89-22_150.txt");
    ax = L1_helper(tmp,maxMld,thresh,testSel,"ad",true,1);
    sgtitle("Low-Level Nitrogen 89-22: L1"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L1/bottle/" + lp + "lln" + tmpT + ".png");
    clear tmp ax;
    
    % Bottle Dissolved Oxygen (88-21)
    tmp = importdata("data/L1/oxy_88-21_150.txt");
    ax = L1_helper(tmp,maxMld,thresh,testSel,"ad",true,1);
    sgtitle("Dissolved Oxygen 88-21: L1"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L1/bottle/" + lp + "boxy" + tmpT + ".png");
    clear tmp ax;
    
    % PProd Light-12 (89-22) (4 significant digits => very good!)
    tmp = importdata("data\L1\l12_89-22_150.txt");
    ax = L1_helper(tmp,maxMld,thresh,testSel,"ad",true,1);
    sgtitle("PProd Light-12 89-22: L1"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L1/bottle/" + lp + "l12" + tmpT + ".png");
    clear tmp ax;

    %%% SPRING
    tmpT = "-ad-02";

    % chla
    tmpx = " Spring";
    tmp = importdata("data/L1/hplcChla_88-21_150.txt");
    ax = L1_helper(tmp,maxMld,thresh,testSel,"ad",true,2);
    sgtitle("L1" + tmpx,"Interpreter","latex");
    exportgraphics(ax,"figures/L1/bottle/" + lp + "chla" + tmpT + ".png");
    clear tmp ax;

    % mvchla (94-21)
    tmp = importdata("data\L1\mvchla_94-21_150.txt");
    ax = L1_helper(tmp,maxMld,thresh,testSel,"ad",true,2);
    sgtitle("HPLC Monovinyl chlorophyll a 94-21: L1"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L1/bottle/" + lp + "mvchla" + tmpT + ".png");
    clear tmp ax;
    
    % dvchla (94-21)
    tmp = importdata("data\L1\dvchla_94-21_150.txt");
    ax = L1_helper(tmp,maxMld,thresh,testSel,"ad",true,2);
    sgtitle("HPLC Divinyl chlorophyll a 94-21: L1"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L1/bottle/" + lp + "dvchla" + tmpT + ".png");
    clear tmp ax;
    
    % chlb (88-21)
    tmp = importdata("data\L1\chlb_88-21_150.txt");
    ax = L1_helper(tmp,maxMld,thresh,testSel,"ad",true,2);
    sgtitle("HPLC chlorophyll b 88-21: L1"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L1/bottle/" + lp + "chlb" + tmpT + ".png");
    clear tmp ax;
    
    % HPLC Chlorophyll C1 + C2 + C3 (88-21)
    tmp = importdata("data\L1\chl123_88-21_150.txt");
    ax = L1_helper(tmp,maxMld,thresh,testSel,"ad",true,2);
    sgtitle("HPLC Chl c1 + c2 + c3: L1"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L1/bottle/" + lp + "chlc123" + tmpT + ".png");
    clear tmp ax;
    
    % HPLC alpha-Carotene (94-21)
    tmp = importdata("data\L1\acar_94-21_150.txt");
    ax = L1_helper(tmp,maxMld,thresh,testSel,"ad",true,2);
    sgtitle("HPLC alpha-Carotene 94-21: L1"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L1/bottle/" + lp + "acar" + tmpT + ".png");
    clear tmp ax;
    
    % HPLC 19' Butanoyloxyfucoxanthin (88-21)
    tmp = importdata("data\L1\but19_88-21_150.txt");
    ax = L1_helper(tmp,maxMld,thresh,testSel,"ad",true,2);
    sgtitle("HPLC 19 Butanoyloxyfucoxanthin 88-21: L1"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L1/bottle/" + lp + "but19" + tmpT + ".png");
    clear tmp ax;
    
    % HPLC 19" Hexanoyloxyfucoxanthin (88-21)
    tmp = importdata("data\L1\hex19_88-21_150.txt");
    ax = L1_helper(tmp,maxMld,thresh,testSel,"ad",true,2);
    sgtitle("HPLC 19' Hexanoyloxyfucoxanthin 88-21: L1"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L1/bottle/" + lp + "hex19" + tmpT + ".png");
    clear tmp ax;
    
    % HPLC Zeaxanthin (88-21)
    tmp = importdata("data\L1\zeaxan_88-21_150.txt");
    ax = L1_helper(tmp,maxMld,thresh,testSel,"ad",true,2);
    sgtitle("HPLC Zeaxanthin 88-21: L1"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L1/bottle/" + lp + "zeax" + tmpT + ".png");
    clear tmp ax;
    
    % Particulate Carbon (89-21)
    tmp = importdata("data\L1\parc_89-21_150.txt");
    ax = L1_helper(tmp,maxMld,thresh,testSel,"ad",true,2);
    sgtitle("L1 Spring","Interpreter","latex");
    exportgraphics(ax,"figures/L1/bottle/" + lp + "pc" + tmpT + ".png");
    clear tmp ax;
    
    % Particulate Nitrogen (89-21)
    tmp = importdata("data\L1\parn_89-21_150.txt");
    ax = L1_helper(tmp,maxMld,thresh,testSel,"ad",true,2);
    sgtitle("Particulate Nitrogen 89-21: L1"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L1/bottle/" + lp + "pn" + tmpT + ".png");
    clear tmp ax;
    
    % Low-Level Phosphorus (88-22)
    tmp = importdata("data\L1\llp_88-22_150.txt");
    ax = L1_helper(tmp,maxMld,thresh,testSel,"ad",true,2);
    sgtitle("Low-Level Phosphorus 88-22: L1"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L1/bottle/" + lp + "llp" + tmpT + ".png");
    clear tmp ax;
    
    % Low-Level Nitrogen (89-22)
    tmp = importdata("data\L1\lln_89-22_150.txt");
    ax = L1_helper(tmp,maxMld,thresh,testSel,"ad",true,2);
    sgtitle("Low-Level Nitrogen 89-22: L1"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L1/bottle/" + lp + "lln" + tmpT + ".png");
    clear tmp ax;
    
    % Bottle Dissolved Oxygen (88-21)
    tmp = importdata("data/L1/oxy_88-21_150.txt");
    ax = L1_helper(tmp,maxMld,thresh,testSel,"ad",true,2);
    sgtitle("Dissolved Oxygen 88-21: L1"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L1/bottle/" + lp + "boxy" + tmpT + ".png");
    clear tmp ax;
    
    % PProd Light-12 (89-22) (4 significant digits => very good!)
    tmp = importdata("data\L1\l12_89-22_150.txt");
    ax = L1_helper(tmp,maxMld,thresh,testSel,"ad",true,2);
    sgtitle("PProd Light-12 89-22: L1"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L1/bottle/" + lp + "l12" + tmpT + ".png");
    clear tmp ax;

    %%% SUMMER
    tmpT = "-ad-03";

    % chla
    tmpx = "Summer";
    tmp = importdata("data/L1/hplcChla_88-21_150.txt");
    ax = L1_helper(tmp,maxMld,thresh,testSel,"ad",true,3);
    sgtitle("L1 " + tmpx,"Interpreter","latex");
    exportgraphics(ax,"figures/L1/bottle/" + lp + "chla" + tmpT + ".png");
    clear tmp ax;

    % mvchla (94-21)
    tmp = importdata("data\L1\mvchla_94-21_150.txt");
    ax = L1_helper(tmp,maxMld,thresh,testSel,"ad",true,3);
    sgtitle("HPLC Monovinyl chlorophyll a 94-21: L1"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L1/bottle/" + lp + "mvchla" + tmpT + ".png");
    clear tmp ax;
    
    % dvchla (94-21)
    tmp = importdata("data\L1\dvchla_94-21_150.txt");
    ax = L1_helper(tmp,maxMld,thresh,testSel,"ad",true,3);
    sgtitle("HPLC Divinyl chlorophyll a 94-21: L1"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L1/bottle/" + lp + "dvchla" + tmpT + ".png");
    clear tmp ax;
    
    % chlb (88-21)
    tmp = importdata("data\L1\chlb_88-21_150.txt");
    ax = L1_helper(tmp,maxMld,thresh,testSel,"ad",true,3);
    sgtitle("HPLC chlorophyll b 88-21: L1"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L1/bottle/" + lp + "chlb" + tmpT + ".png");
    clear tmp ax;
    
    % HPLC Chlorophyll C1 + C2 + C3 (88-21)
    tmp = importdata("data\L1\chl123_88-21_150.txt");
    ax = L1_helper(tmp,maxMld,thresh,testSel,"ad",true,3);
    sgtitle("HPLC Chl c1 + c2 + c3: L1"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L1/bottle/" + lp + "chlc123" + tmpT + ".png");
    clear tmp ax;
    
    % HPLC alpha-Carotene (94-21)
    tmp = importdata("data\L1\acar_94-21_150.txt");
    ax = L1_helper(tmp,maxMld,thresh,testSel,"ad",true,3);
    sgtitle("HPLC alpha-Carotene 94-21: L1"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L1/bottle/" + lp + "acar" + tmpT + ".png");
    clear tmp ax;
    
    % HPLC 19' Butanoyloxyfucoxanthin (88-21)
    tmp = importdata("data\L1\but19_88-21_150.txt");
    ax = L1_helper(tmp,maxMld,thresh,testSel,"ad",true,3);
    sgtitle("HPLC 19 Butanoyloxyfucoxanthin 88-21: L1"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L1/bottle/" + lp + "but19" + tmpT + ".png");
    clear tmp ax;
    
    % HPLC 19" Hexanoyloxyfucoxanthin (88-21)
    tmp = importdata("data\L1\hex19_88-21_150.txt");
    ax = L1_helper(tmp,maxMld,thresh,testSel,"ad",true,3);
    sgtitle("HPLC 19' Hexanoyloxyfucoxanthin 88-21: L1"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L1/bottle/" + lp + "hex19" + tmpT + ".png");
    clear tmp ax;
    
    % HPLC Zeaxanthin (88-21)
    tmp = importdata("data\L1\zeaxan_88-21_150.txt");
    ax = L1_helper(tmp,maxMld,thresh,testSel,"ad",true,3);
    sgtitle("HPLC Zeaxanthin 88-21: L1"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L1/bottle/" + lp + "zeax" + tmpT + ".png");
    clear tmp ax;
    
    % Particulate Carbon (89-21)
    tmp = importdata("data\L1\parc_89-21_150.txt");
    ax = L1_helper(tmp,maxMld,thresh,testSel,"ad",true,3);
    sgtitle("L1 Summer","Interpreter","latex");
    exportgraphics(ax,"figures/L1/bottle/" + lp + "pc" + tmpT + ".png");
    clear tmp ax;
    
    % Particulate Nitrogen (89-21)
    tmp = importdata("data\L1\parn_89-21_150.txt");
    ax = L1_helper(tmp,maxMld,thresh,testSel,"ad",true,3);
    sgtitle("Particulate Nitrogen 89-21: L1"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L1/bottle/" + lp + "pn" + tmpT + ".png");
    clear tmp ax;
    
    % Low-Level Phosphorus (88-22)
    tmp = importdata("data\L1\llp_88-22_150.txt");
    ax = L1_helper(tmp,maxMld,thresh,testSel,"ad",true,3);
    sgtitle("Low-Level Phosphorus 88-22: L1"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L1/bottle/" + lp + "llp" + tmpT + ".png");
    clear tmp ax;
    
    % Low-Level Nitrogen (89-22)
    tmp = importdata("data\L1\lln_89-22_150.txt");
    ax = L1_helper(tmp,maxMld,thresh,testSel,"ad",true,3);
    sgtitle("Low-Level Nitrogen 89-22: L1"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L1/bottle/" + lp + "lln" + tmpT + ".png");
    clear tmp ax;
    
    % Bottle Dissolved Oxygen (88-21)
    tmp = importdata("data/L1/oxy_88-21_150.txt");
    ax = L1_helper(tmp,maxMld,thresh,testSel,"ad",true,3);
    sgtitle("Dissolved Oxygen 88-21: L1"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L1/bottle/" + lp + "boxy" + tmpT + ".png");
    clear tmp ax;
    
    % PProd Light-12 (89-22) (4 significant digits => very good!)
    tmp = importdata("data\L1\l12_89-22_150.txt");
    ax = L1_helper(tmp,maxMld,thresh,testSel,"ad",true,3);
    sgtitle("PProd Light-12 89-22: L1"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L1/bottle/" + lp + "l12" + tmpT + ".png");
    clear tmp ax;

    %%% AUTUMN
    tmpT = "-ad-04";

    % chla
    tmpx = " Autumn";
    tmp = importdata("data/L1/hplcChla_88-21_150.txt");
    ax = L1_helper(tmp,maxMld,thresh,testSel,"ad",true,4);
    sgtitle("L1" + tmpx,"Interpreter","latex");
    exportgraphics(ax,"figures/L1/bottle/" + lp + "chla" + tmpT + ".png");
    clear tmp ax;

    % mvchla (94-21)
    tmp = importdata("data\L1\mvchla_94-21_150.txt");
    ax = L1_helper(tmp,maxMld,thresh,testSel,"ad",true,4);
    sgtitle("HPLC Monovinyl chlorophyll a 94-21: L1"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L1/bottle/" + lp + "mvchla" + tmpT + ".png");
    clear tmp ax;
    
    % dvchla (94-21)
    tmp = importdata("data\L1\dvchla_94-21_150.txt");
    ax = L1_helper(tmp,maxMld,thresh,testSel,"ad",true,4);
    sgtitle("HPLC Divinyl chlorophyll a 94-21: L1"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L1/bottle/" + lp + "dvchla" + tmpT + ".png");
    clear tmp ax;
    
    % chlb (88-21)
    tmp = importdata("data\L1\chlb_88-21_150.txt");
    ax = L1_helper(tmp,maxMld,thresh,testSel,"ad",true,4);
    sgtitle("HPLC chlorophyll b 88-21: L1"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L1/bottle/" + lp + "chlb" + tmpT + ".png");
    clear tmp ax;
    
    % HPLC Chlorophyll C1 + C2 + C3 (88-21)
    tmp = importdata("data\L1\chl123_88-21_150.txt");
    ax = L1_helper(tmp,maxMld,thresh,testSel,"ad",true,4);
    sgtitle("HPLC Chl c1 + c2 + c3: L1"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L1/bottle/" + lp + "chlc123" + tmpT + ".png");
    clear tmp ax;
    
    % HPLC alpha-Carotene (94-21)
    tmp = importdata("data\L1\acar_94-21_150.txt");
    ax = L1_helper(tmp,maxMld,thresh,testSel,"ad",true,4);
    sgtitle("HPLC alpha-Carotene 94-21: L1"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L1/bottle/" + lp + "acar" + tmpT + ".png");
    clear tmp ax;
    
    % HPLC 19' Butanoyloxyfucoxanthin (88-21)
    tmp = importdata("data\L1\but19_88-21_150.txt");
    ax = L1_helper(tmp,maxMld,thresh,testSel,"ad",true,4);
    sgtitle("HPLC 19 Butanoyloxyfucoxanthin 88-21: L1"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L1/bottle/" + lp + "but19" + tmpT + ".png");
    clear tmp ax;
    
    % HPLC 19" Hexanoyloxyfucoxanthin (88-21)
    tmp = importdata("data\L1\hex19_88-21_150.txt");
    ax = L1_helper(tmp,maxMld,thresh,testSel,"ad",true,4);
    sgtitle("HPLC 19' Hexanoyloxyfucoxanthin 88-21: L1"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L1/bottle/" + lp + "hex19" + tmpT + ".png");
    clear tmp ax;
    
    % HPLC Zeaxanthin (88-21)
    tmp = importdata("data\L1\zeaxan_88-21_150.txt");
    ax = L1_helper(tmp,maxMld,thresh,testSel,"ad",true,4);
    sgtitle("HPLC Zeaxanthin 88-21: L1"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L1/bottle/" + lp + "zeax" + tmpT + ".png");
    clear tmp ax;
    
    % Particulate Carbon (89-21)
    tmp = importdata("data\L1\parc_89-21_150.txt");
    ax = L1_helper(tmp,maxMld,thresh,testSel,"ad",true,4);
    sgtitle("L1 Autumn","Interpreter","latex");
    exportgraphics(ax,"figures/L1/bottle/" + lp + "pc" + tmpT + ".png");
    clear tmp ax;
    
    % Particulate Nitrogen (89-21)
    tmp = importdata("data\L1\parn_89-21_150.txt");
    ax = L1_helper(tmp,maxMld,thresh,testSel,"ad",true,4);
    sgtitle("Particulate Nitrogen 89-21: L1"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L1/bottle/" + lp + "pn" + tmpT + ".png");
    clear tmp ax;
    
    % Low-Level Phosphorus (88-22)
    tmp = importdata("data\L1\llp_88-22_150.txt");
    ax = L1_helper(tmp,maxMld,thresh,testSel,"ad",true,4);
    sgtitle("Low-Level Phosphorus 88-22: L1"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L1/bottle/" + lp + "llp" + tmpT + ".png");
    clear tmp ax;
    
    % Low-Level Nitrogen (89-22)
    tmp = importdata("data\L1\lln_89-22_150.txt");
    ax = L1_helper(tmp,maxMld,thresh,testSel,"ad",true,4);
    sgtitle("Low-Level Nitrogen 89-22: L1"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L1/bottle/" + lp + "lln" + tmpT + ".png");
    clear tmp ax;
    
    % Bottle Dissolved Oxygen (88-21)
    tmp = importdata("data/L1/oxy_88-21_150.txt");
    ax = L1_helper(tmp,maxMld,thresh,testSel,"ad",true,4);
    sgtitle("Dissolved Oxygen 88-21: L1"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L1/bottle/" + lp + "boxy" + tmpT + ".png");
    clear tmp ax;
    
    % PProd Light-12 (89-22) (4 significant digits => very good!)
    tmp = importdata("data\L1\l12_89-22_150.txt");
    ax = L1_helper(tmp,maxMld,thresh,testSel,"ad",true,4);
    sgtitle("PProd Light-12 89-22: L1"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L1/bottle/" + lp + "l12" + tmpT + ".png");
    clear tmp ax;


end

%% Parameters not used for current analysis

% % Dissolved Organic Nitrogen (88-17)
% tmp = importdata("data\L1\don_88-17_150.txt");
% [ax,p,ks,obs,Sk,Ku] = L1_helper(tmp,maxMld);
% sgtitle("DON 88-17: L1");
% exportgraphics(ax,"figures/L1/bottle/" + lp + "notUsed/don" + tmpT + ".png");
% save("output\L1\don.mat","p","ks","obs","Sk","Ku");
% clearvars -except tmpT maxMld lp;

% % Dissolved Organic Carbon (DOC): 93-17
% tmp = importdata("data\L1\doc_93-17_150.txt");
% [ax,p,ks,obs,Sk,Ku] = L1_helper(tmp,maxMld);
% sgtitle("DOC 93-17: L1");
% exportgraphics(ax,"figures/L1/bottle/" + lp + "notUsed/doc" + tmpT + ".png");
% save("output\L1\doc.mat","p","ks","obs","Sk","Ku");
% clearvars -except tmpT maxMld lp;

% % Total Dissolved Phosphorus (TDP): 88-01
% tmp = importdata("data\L1\tdp_88-01_150.txt");
% [ax,p,ks,obs,Sk,Ku] = L1_helper(tmp,maxMld);
% sgtitle("Total Dissolved Phosphorus 88-01: L1");
% exportgraphics(ax,"figures/L1/bottle/" + lp + "notUsed/tdp" + tmpT + ".png");
% save("output\L1\tdp.mat","p","ks","obs","Sk","Ku");
% clearvars -except tmpT maxMld lp;

% % Total Dissolved Nitrogen: 88-17
% tmp = importdata("data\L1\tdn_88-17_150.txt");
% [ax,p,ks,obs,Sk,Ku] = L1_helper(tmp,maxMld);
% sgtitle("Total Dissolved Nitrogen 88-17: L1");
% exportgraphics(ax,"figures/L1/bottle/" + lp + "notUsed/tdn" + tmpT + ".png");
% save("output\L1\tdn.mat","p","ks","obs","Sk","Ku");
% clearvars -except tmpT maxMld lp;

% % Particulate Phosphorus: 11-21
% tmp = importdata("data\L1\parp_11-21_150.txt");
% [ax,p,ks,obs,Sk,Ku] = L1_helper(tmp,maxMld);
% sgtitle("Particulate Phosphorus 11-21: L1");
% exportgraphics(ax,"figures/L1/bottle/" + lp + "notUsed/pp" + tmpT + ".png");
% save("output\L1\pp.mat","p","ks","obs","Sk","Ku");
% clearvars -except tmpT maxMld lp;

% % Fluorometric Chlorophyll a (89-22)
% tmp = importdata("data\L1\chlFlu_88-22_150.txt");
% [ax,p,ks,obs,Sk,Ku] = L1_helper(tmp,maxMld);
% sgtitle("Fluorometric Chlorophyll a 89-22: L1");
% exportgraphics(ax,"figures/L1/bottle/" + lp + "notUsed/chlaFluo" + tmpT + ".png");
% save("output\L1\chlaFluo.mat","p","ks","obs","Sk","Ku");
% clearvars -except tmpT maxMld lp;

% % Phaeopigments (88-22)
% tmp = importdata("data\L1\pheo_88-22_150.txt");
% [ax,p,ks,obs,Sk,Ku] = L1_helper(tmp,maxMld);
% sgtitle("Phaeopigments 88-22: L1");
% exportgraphics(ax,"figures/L1/bottle/" + lp + "notUsed/phaeo" + tmpT + ".png");
% save("output\L1\phaeo.mat","p","ks","obs","Sk","Ku");
% clearvars -except tmpT maxMld lp;

% % HPLC Chlorophyll C3 (88-21)
% tmp = importdata("data\L1\chl3_88-21_150.txt");
% [ax,p,ks,obs,Sk,Ku] = L1_helper(tmp,maxMld);
% sgtitle("HPLC Chl c3: L1");
% exportgraphics(ax,"figures/L1/bottle/" + lp + "notUsed/chl3" + tmpT + ".png");
% save("output\L1\chl3.mat","p","ks","obs","Sk","Ku");
% clearvars -except tmpT maxMld lp;

% % HPLC Chlorophyll C1 + C2: 88-21
% tmp = importdata("data\L1\chl12_88-21_150.txt");
% [ax,p,ks,obs,Sk,Ku] = L1_helper(tmp,maxMld); 
% sgtitle("HPLC Chl c1 + c2: L1");
% exportgraphics(ax,"figures/L1/bottle/" + lp + "notUsed/chl12" + tmpT + ".png");
% save("output\L1\chl12.mat","p","ks","obs","Sk","Ku");
% clearvars -except tmpT maxMld lp;

% % Peridinin (88-21)
% tmp = importdata("data\L1\per_88-21_150.txt");
% [ax,p,ks,obs,Sk,Ku] = L1_helper(tmp,maxMld);
% sgtitle("Peridinin: L1");
% exportgraphics(ax,"figures/L1/bottle/" + lp + "notUsed/perid" + tmpT + ".png");
% save("output\L1\perid.mat","p","ks","obs","Sk","Ku");
% clearvars -except tmpT maxMld lp;

% % HPLC Fucoxanthin (88-21)
% tmp = importdata("data\L1\fuco_88-21_150.txt");
% [ax,p,ks,obs,Sk,Ku] = L1_helper(tmp,maxMld);
% sgtitle("Fucoxanthin 88-21: L1");
% exportgraphics(ax,"figures/L1/bottle/" + lp + "notUsed/fuco" + tmpT + ".png"); clear ax;
% save("output\L1\fuco.mat","p","ks","obs","Sk","Ku");
% clearvars -except tmpT maxMld lp;

% % HPLC Diadinoxanthin (88-21)
% tmp = importdata("data\L1\diadino_88-21_150.txt");
% [ax,p,ks,obs,Sk,Ku] = L1_helper(tmp,maxMld);
% sgtitle("HPLC Diadinoxanthin 88-21: L1");
% exportgraphics(ax,"figures/L1/bottle/" + lp + "notUsed/diadino" + tmpT + ".png");
% save("output\L1\diadino.mat","p","ks","obs","Sk","Ku");
% clearvars -except tmpT maxMld lp;

% % HPLC beta-Carotene (94-21)
% tmp = importdata("data\L1\bcar_94-21_150.txt");
% [ax,p,ks,obs,Sk,Ku] = L1_helper(tmp,maxMld);
% sgtitle("HPLC beta-Carotene 94-21: L1");
% exportgraphics(ax,"figures/L1/bottle/" + lp + "notUsed/bcar" + tmpT + ".png");
% save("output\L1\bcar.mat","p","ks","obs","Sk","Ku");
% clearvars -except tmpT maxMld lp;

% % HPLC Carotenes (88-21)
% tmp = importdata("data\L1\caroten_88-21_150.txt");
% [ax,p,ks,obs,Sk,Ku] = L1_helper(tmp,maxMld);
% sgtitle("HPLC Carotenes 88-21: L1");
% exportgraphics(ax,"figures/L1/bottle/" + lp + "notUsed/caroten" + tmpT + ".png"); clear ax;
% save("output\L1\caroten.mat","p","ks","obs","Sk","Ku");

% % HPLC Chlorophyllide a (94-21)
% tmp = importdata("data\L1\chlda_94-21_150.txt");
% [ax,p,ks,obs,Sk,Ku] = L1_helper(tmp,maxMld);
% sgtitle("HPLC Chlorophyllide a 94-21: L1");
% exportgraphics(ax,"figures/L1/bottle/" + lp + "notUsed/chlda" + tmpT + ".png");
% save("output\L1\chlda.mat","p","ks","obs","Sk","Ku");
% clearvars -except tmpT maxMld lp;

% % HPLC Lutein (94-21)
% tmp = importdata("data\L1\lutein_94-21_150.txt");
% [ax,p,ks,obs,Sk,Ku] = L1_helper(tmp,maxMld);
% sgtitle("HPLC Lutein 94-21: L1");
% exportgraphics(ax,"figures/L1/bottle/" + lp + "notUsed/lutein" + tmpT + ".png");
% save("output\L1\lutein.mat","p","ks","obs","Sk","Ku");
% clearvars -except tmpT maxMld lp;

% % Phycoerythrin 0.4u fraction (00-08)
% tmp = importdata("data\L1\pe4_00-08_150.txt");
% [ax,p,ks,obs,Sk,Ku] = L1_helper(tmp,maxMld);
% sgtitle("Phycoerythrin 0.4u fraction 00-08: L1");
% exportgraphics(ax,"figures/L1/bottle/" + lp + "notUsed/pe4" + tmpT + ".png");
% save("output\L1\pe4.mat","p","ks","obs","Sk","Ku");
% clearvars -except tmpT maxMld lp;

% % Phycoerythrin 5u fraction (00-08) (Three significant digits => good!)
% tmp = importdata("data\L1\pe5_00-08_150.txt");
% [ax,p,ks,obs,Sk,Ku] = L1_helper(tmp,maxMld);
% sgtitle("Phycoerythrin 5u fraction 00-08: L1");
% exportgraphics(ax,"figures/L1/bottle/" + lp + "notUsed/pe5" + tmpT + ".png");
% save("output\L1\pe5.mat","p","ks","obs","Sk","Ku");
% clearvars -except tmpT maxMld lp;

% % Phycoerythrin 10u fraction (00-08) (4 significant digits => very good!)
% tmp = importdata("data\L1\pe10_00-08_150.txt");
% [ax,p,ks,obs,Sk,Ku] = L1_helper(tmp,maxMld);
% sgtitle("Phycoerythrin 10u fraction 00-08: L1");
% exportgraphics(ax,"figures/L1/bottle/" + lp + "notUsed/pe10" + tmpT + ".png");
% save("output\L1\pe10.mat","p","ks","obs","Sk","Ku");
% clearvars -except tmpT maxMld lp;

% % Heterotrophic Bacteria (05-21) (4 significant digits => very good!)
% tmp = importdata("data\L1\hbact_05-21_150.txt");
% [ax,p,ks,obs,Sk,Ku] = L1_helper(tmp,maxMld);
% sgtitle("Heterotrophic Bacteria 05-21: L1");
% exportgraphics(ax,"figures/L1/bottle/" + lp + "notUsed/hbact" + tmpT + ".png");
% save("output\L1\hbact.mat","p","ks","obs","Sk","Ku");
% clearvars -except tmpT maxMld lp;

% % Prochlorococcus (05-21) (4 significant digits => very good!)
% tmp = importdata("data\L1\pbact_05-21_150.txt");
% [ax,p,ks,obs,Sk,Ku] = L1_helper(tmp,maxMld);
% sgtitle("Prochlorococcus 05-21: L1");
% exportgraphics(ax,"figures/L1/pbact.png"); clear ax;
% save("output\L1\pbact.mat","p","ks","obs","Sk","Ku");
% clearvars -except tmpT maxMld lp;

% % Synechococcus (05-21) (Two significant digits => may be unreliable!)
% tmp = importdata("data\L1\sbact_05-21_150.txt");
% [ax,p,ks,obs,Sk,Ku] = L1_helper(tmp,maxMld);
% sgtitle("Synechococcus 05-21: L1");
% exportgraphics(ax,"figures/L1/bottle/" + lp + "notUsed/sbact" + tmpT + ".png");
% save("output\L1\sbact.mat","p","ks","obs","Sk","Ku");
% clearvars -except tmpT maxMld lp;

% % PicoEukaryotes: 05-21 (2 significant digits => may be unreliable!)
% tmp = importdata("data\L1\ebact_05-21_150.txt");
% [ax,p,ks,obs,Sk,Ku] = L1_helper(tmp,maxMld);
% sgtitle("Picoeukaryotes 05-21: L1");
% exportgraphics(ax,"figures/L1/bottle/" + lp + "notUsed/ebact" + tmpT + ".png");
% save("output\L1\ebact.mat","p","ks","obs","Sk","Ku");
% clearvars -except tmpT maxMld lp;

% % ATP (88-22) (Four significant digits => very good!)
% tmp = importdata("data\L1\atp_88-22_150.txt");
% [ax,p,ks,obs,Sk,Ku] = L1_helper(tmp,maxMld);
% sgtitle("ATP 88-22: L1");
% exportgraphics(ax,"figures/L1/bottle/" + lp + "notUsed/atp" + tmpT + ".png");
% save("output\L1\atp.mat","p","ks","obs","Sk","Ku");
% clearvars -except tmpT maxMld lp;

% % PProd Dark-12 (89-00) (Three significant digits => good!)
% tmp = importdata("data\L1\d12_89-00_150.txt");
% [ax,p,ks,obs,Sk,Ku] = L1_helper(tmp,maxMld);
% sgtitle("PProd Dark-12 89-00: L1");
% exportgraphics(ax,"figures/L1/bottle/" + lp + "notUsed/d12" + tmpT + ".png");
% save("output\L1\d12.mat","p","ks","obs","Sk","Ku");
% clearvars -except tmpT maxMld lp;

% Dissolved Organic Phosphate 
% tmp = importdata("data\L1\pho_88-21_150.txt");
% [ax,p,ks,obs,Sk,Ku] = L1_helper(tmp,maxMld);
% sgtitle("Phosphate 88-21: L1");
% exportgraphics(ax,"figures/L1/bottle/phos" + tmpT + ".png");
% save("output\L1\phos.mat","p","ks","obs","Sk","Ku");
% clearvars -except tmpT maxMld lp;