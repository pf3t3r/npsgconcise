% Script to output L2 bottle results for the statistical analysis.

clear; clc; close all;
addpath("baroneRoutines\"); addpath("func\");
set(groot, "defaultFigureUnits", "centimeters", "defaultFigurePosition", [3 3 15 15]);

% Load maximum mixed-layer depth 'MLD' and cruise-averaged deep chlorophyll
% maximum depth 'DCM'.
mld = load("mldVals.mat").maxMld; % single maximum per cruise
dcm = load("output/dcm.mat").dcm; % pDcm + sigmaDcm (all casts, crn 1-329)

% Possible test cases.
principleAnalysis = false;  % main analysis
seasonalAnalysisKs = false;   % seasonality of statistics
seasonalAnalysisAd = true;   % seasonality of statistics
testSel = 2; % 2 = norm + logn; 4 = norm + logn + weib + gamm

% TEMPLATE (XX-YY)
% 1. Load data
% 2. Extract data beneath ML, centre around DCM
% 3. Calculate K-S (or A-D) p-value, Vuong LLR, skewness + kurtosis
% 4. Plot results

%% L2 Histograms
tmp = importdata("data/L2/hplcChla_88-21_200.txt");
% [~,~,pSubml,~,~,~,~,~,~,~,~,~,~,~,cSubml] = L2_helper(tmp,mld,dcm,30,2,"ad");
% pB = round(pSubml,-1);   % bin the pressure

[~,chla_L2,~,~,~,~,~,~,~,~,~,~,~,~,~,p_L2] = L2_helper(tmp,mld,dcm,30,2,"ad");

%% hist

figure
histfit(chla_L2(p_L2==60),10,"lognormal");

%% Seasonal Analysis: K-S
% thresh = 30;
if seasonalAnalysisKs == true
    
    thresh = 50;
    %%% WINTER
    tmpT = "-01";

    % Chlorophyll a (88-21)
    tmp = importdata("data/L2/hplcChla_88-21_200.txt");
    ax = L2_helper(tmp,mld,dcm,thresh,testSel,"ks",[-60 60],[5 17],1);
    sgtitle("[Chl a] 88-21: L2"+tmpT);
    exportgraphics(ax,"figures/L2/bottle/log/chla" + tmpT + ".png");
    clear tmp ax;

    % HPLC Monovinyl chlorophyll a (94-21)
    tmp = importdata("data\L2\mvchla_94-21_200.txt");
    ax = L2_helper(tmp,mld,dcm,thresh,testSel,"ks",[-60 60],[5 17],1);
    sgtitle("HPLC Monovinyl chlorophyll a (94-21): L2"+tmpT);
    exportgraphics(ax,"figures/L2/bottle/log/mvchla" + tmpT + ".png");
    clear tmp ax;
    
    % HPLC Divinyl chlorophyll a (94-21)
    tmp = importdata("data\L2\dvchla_94-21_200.txt");
    ax = L2_helper(tmp,mld,dcm,thresh,testSel,"ks",[-60 60],[5 17],1);
    sgtitle("HPLC Divinyl chlorophyll a (94-21): L2"+tmpT);
    exportgraphics(ax,"figures/L2/bottle/log/dvchla" + tmpT + ".png");
    clear tmp ax;
    
    % HPLC Chlorophyll b (88-21)
    tmp = importdata("data\L2\chlb_88-21_200.txt");
    ax = L2_helper(tmp,mld,dcm,thresh,testSel,"ks",[-60 60],[5 17],1);
    sgtitle("HPLC Chlorophyll b (88-21): L2"+tmpT);
    exportgraphics(ax,"figures/L2/bottle/log/chlb" + tmpT + ".png");
    clear tmp ax;
    
    % HPLC Chlorophyll C1 + C2 + C3 (88-21)
    tmp = importdata("data\L2\chl123_88-21_200.txt");
    ax = L2_helper(tmp,mld,dcm,thresh,testSel,"ks",[-60 60],[5 17],1);
    sgtitle("Chlorophyll C1 + C2 + C3 (88-21): L2"+tmpT);
    exportgraphics(ax,"figures/L2/bottle/log/chl123" + tmpT + ".png");
    clear tmp ax;
    
    % HPLC alpha-Carotene (94-21)
    tmp = importdata("data\L2\acar_94-21_200.txt");
    ax = L2_helper(tmp,mld,dcm,thresh,testSel,"ks",[-60 60],[5 17],1);
    sgtitle("HPLC alpha-Carotene (94-21): L2"+tmpT);
    exportgraphics(ax,"figures/L2/bottle/log/acar" + tmpT + ".png");
    clear tmp ax;
    
    % But19 (88-21)
    tmp = importdata("data\L2\but19_88-21_200.txt");
    ax = L2_helper(tmp,mld,dcm,thresh,testSel,"ks",[-60 60],[5 17],1);
    sgtitle("But19 (88-21): L2" + tmpT);
    exportgraphics(ax,"figures/L2/bottle/log/but19" + tmpT + ".png");
    clear tmp ax;
    
    % HPLC 19 Hexanoyloxyfucoxanthin (88-21)
    tmp = importdata("data\L2\hex19_88-21_200.txt");
    ax = L2_helper(tmp,mld,dcm,thresh,testSel,"ks",[-60 60],[5 17],1);
    sgtitle("HPLC 19' Hexanoyloxyfucoxanthin (88-21): L2" + tmpT);
    exportgraphics(ax,"figures/L2/bottle/log/hex19" + tmpT + ".png");
    clear tmp ax;
    
    % HPLC Zeaxanthin (88-21)
    tmp = importdata("data\L2\zeaxan_88-21_200.txt");
    ax = L2_helper(tmp,mld,dcm,thresh,testSel,"ks",[-60 60],[5 17],1);
    sgtitle("HPLC Zeaxanthin (88-21): L2" + tmpT);
    exportgraphics(ax,"figures/L2/bottle/log/zeaxan" + tmpT + ".png");
    clear tmp ax;
    
    % Particulate Carbon (89-21)
    tmp = importdata("data\L2\parc_89-21_200.txt");
    ax = L2_helper(tmp,mld,dcm,thresh,testSel,"ks",[-60 60],[2 14],1);
    sgtitle("Particulate Carbon 89-21: L2" + tmpT);
    exportgraphics(ax,"figures/L2/bottle/log/pc" + tmpT + ".png");
    clear tmp ax;
    
    % Particulate Nitrogen (89-21)
    tmp = importdata("data\L2\parn_89-21_200.txt");
    ax = L2_helper(tmp,mld,dcm,thresh,testSel,"ks",[-60 60],[2 14],1);
    sgtitle("Particulate Nitrogen 89-21: L2" + tmpT);
    exportgraphics(ax,"figures/L2/bottle/log/pn" + tmpT + ".png");
    clear tmp ax;
    
    % Nitrate + Nitrite (88-21)
    tmp = importdata("data/L2/nit2_88-21_200.txt");
    ax = L2_helper(tmp,mld,dcm,thresh,testSel,"ks",[-60 60],[6 18],1);
    sgtitle("Nitrate + Nitrite 88-21: L2" + tmpT);
    exportgraphics(ax,"figures/L2/bottle/log/nit2" + tmpT + ".png");
    clear tmp ax;
    
    % Phosphate (88-21)
    tmp = importdata("data/L2/pho_88-21_200.txt");
    ax = L2_helper(tmp,mld,dcm,thresh,testSel,"ks",[-60 60],[6 18],1);
    sgtitle("Phosphate 88-21: L2" + tmpT);
    exportgraphics(ax,"figures/L2/bottle/log/phos" + tmpT + ".png");
    clear tmp ax;
    
    % Oxygen (88-21)
    tmp = importdata("data/L2/oxy_88-21_200.txt");
    ax = L2_helper(tmp,mld,dcm,thresh,testSel,"ks",[-60 60],[5 17],1);
    sgtitle("Dissolved Oxygen 88-21: L2" + tmpT);
    exportgraphics(ax,"figures/L2/bottle/log/boxy" + tmpT + ".png");
    clear tmp ax;
    
    % PProd Light-12 (89-22)
    tmp = importdata("data\L2\l12_89-22_200.txt");
    ax = L2_helper(tmp,mld,dcm,thresh,testSel,"ks",[-60 60],[3 15],1);
    sgtitle("PProd Light-12 (89-22): L2" + tmpT);
    exportgraphics(ax,"figures/L2/bottle/log/l12" + tmpT + ".png");
    clear tmp ax;


    %%% SPRING
    tmpT = "-02";

    % chla
    tmp = importdata("data/L2/hplcChla_88-21_200.txt");
    ax = L2_helper(tmp,mld,dcm,thresh,testSel,"ks",[-60 60],[7 19],2);
    sgtitle("[Chl a] 88-21: L2"+tmpT);
    exportgraphics(ax,"figures/L2/bottle/log/chla" + tmpT + ".png");
    clear tmp ax;

    % HPLC Monovinyl chlorophyll a (94-21)
    tmp = importdata("data\L2\mvchla_94-21_200.txt");
    ax = L2_helper(tmp,mld,dcm,thresh,testSel,"ks",[-60 60],[7 19],2);
    sgtitle("HPLC Monovinyl chlorophyll a (94-21): L2"+tmpT);
    exportgraphics(ax,"figures/L2/bottle/log/mvchla" + tmpT + ".png");
    clear tmp ax;
    
    % HPLC Divinyl chlorophyll a (94-21)
    tmp = importdata("data\L2\dvchla_94-21_200.txt");
    ax = L2_helper(tmp,mld,dcm,thresh,testSel,"ks",[-60 60],[7 19],2);
    sgtitle("HPLC Divinyl chlorophyll a (94-21): L2"+tmpT);
    exportgraphics(ax,"figures/L2/bottle/log/dvchla" + tmpT + ".png");
    clear tmp ax;
    
    % HPLC Chlorophyll b (88-21)
    tmp = importdata("data\L2\chlb_88-21_200.txt");
    ax = L2_helper(tmp,mld,dcm,thresh,testSel,"ks",[-60 60],[7 19],2);
    sgtitle("HPLC Chlorophyll b (88-21): L2"+tmpT);
    exportgraphics(ax,"figures/L2/bottle/log/chlb" + tmpT + ".png");
    clear tmp ax;
    
    % HPLC Chlorophyll C1 + C2 + C3 (88-21)
    tmp = importdata("data\L2\chl123_88-21_200.txt");
    ax = L2_helper(tmp,mld,dcm,thresh,testSel,"ks",[-60 60],[7 19],2);
    sgtitle("Chlorophyll C1 + C2 + C3 (88-21): L2"+tmpT);
    exportgraphics(ax,"figures/L2/bottle/log/chl123" + tmpT + ".png");
    clear tmp ax;
    
    % HPLC alpha-Carotene (94-21)
    tmp = importdata("data\L2\acar_94-21_200.txt");
    ax = L2_helper(tmp,mld,dcm,thresh,testSel,"ks",[-60 60],[7 19],2);
    sgtitle("HPLC alpha-Carotene (94-21): L2"+tmpT);
    exportgraphics(ax,"figures/L2/bottle/log/acar" + tmpT + ".png");
    clear tmp ax;
    
    % But19 (88-21)
    tmp = importdata("data\L2\but19_88-21_200.txt");
    ax = L2_helper(tmp,mld,dcm,thresh,testSel,"ks",[-60 60],[7 19],2);
    sgtitle("But19 (88-21): L2" + tmpT);
    exportgraphics(ax,"figures/L2/bottle/log/but19" + tmpT + ".png");
    clear tmp ax;
    
    % HPLC 19 Hexanoyloxyfucoxanthin (88-21)
    tmp = importdata("data\L2\hex19_88-21_200.txt");
    ax = L2_helper(tmp,mld,dcm,thresh,testSel,"ks",[-60 60],[7 19],2);
    sgtitle("HPLC 19' Hexanoyloxyfucoxanthin (88-21): L2" + tmpT);
    exportgraphics(ax,"figures/L2/bottle/log/hex19" + tmpT + ".png");
    clear tmp ax;
    
    % HPLC Zeaxanthin (88-21)
    tmp = importdata("data\L2\zeaxan_88-21_200.txt");
    ax = L2_helper(tmp,mld,dcm,thresh,testSel,"ks",[-60 60],[7 19],2);
    sgtitle("HPLC Zeaxanthin (88-21): L2" + tmpT);
    exportgraphics(ax,"figures/L2/bottle/log/zeaxan" + tmpT + ".png");
    clear tmp ax;
    
    % Particulate Carbon (89-21)
    tmp = importdata("data\L2\parc_89-21_200.txt");
    ax = L2_helper(tmp,mld,dcm,thresh,testSel,"ks",[-60 60],[8 20],2);
    sgtitle("Particulate Carbon 89-21: L2" + tmpT);
    exportgraphics(ax,"figures/L2/bottle/log/pc" + tmpT + ".png");
    clear tmp ax;
    
    % Particulate Nitrogen (89-21)
    tmp = importdata("data\L2\parn_89-21_200.txt");
    ax = L2_helper(tmp,mld,dcm,thresh,testSel,"ks",[-60 60],[8 20],2);
    sgtitle("Particulate Nitrogen 89-21: L2" + tmpT);
    exportgraphics(ax,"figures/L2/bottle/log/pn" + tmpT + ".png");
    clear tmp ax;
    
    % Nitrate + Nitrite (88-21)
    tmp = importdata("data/L2/nit2_88-21_200.txt");
    ax = L2_helper(tmp,mld,dcm,thresh,testSel,"ks",[-60 60],[6 18],2);
    sgtitle("Nitrate + Nitrite 88-21: L2" + tmpT);
    exportgraphics(ax,"figures/L2/bottle/log/nit2" + tmpT + ".png");
    clear tmp ax;
    
    % Phosphate (88-21)
    tmp = importdata("data/L2/pho_88-21_200.txt");
    ax = L2_helper(tmp,mld,dcm,thresh,testSel,"ks",[-60 60],[6 18],2);
    sgtitle("Phosphate 88-21: L2" + tmpT);
    exportgraphics(ax,"figures/L2/bottle/log/phos" + tmpT + ".png");
    clear tmp ax;
    
    % Oxygen (88-21)
    tmp = importdata("data/L2/oxy_88-21_200.txt");
    ax = L2_helper(tmp,mld,dcm,thresh,testSel,"ks",[-60 60],[6 18],2);
    sgtitle("Dissolved Oxygen 88-21: L2" + tmpT);
    exportgraphics(ax,"figures/L2/bottle/log/boxy" + tmpT + ".png");
    clear tmp ax;
    
    % PProd Light-12 (89-22)
    tmp = importdata("data\L2\l12_89-22_200.txt");
    ax = L2_helper(tmp,mld,dcm,10,testSel,"ks",[-60 60],[7 19],2);
    sgtitle("PProd Light-12 (89-22): L2" + tmpT);
    exportgraphics(ax,"figures/L2/bottle/log/l12" + tmpT + ".png");
    clear tmp ax;


    %%% SUMMER
    tmpT = "-03";

    % chla
    tmp = importdata("data/L2/hplcChla_88-21_200.txt");
    ax = L2_helper(tmp,mld,dcm,15,testSel,"ks",[-80 80],[3 19],3);
    sgtitle("[Chl a] 88-21: L2"+tmpT);
    exportgraphics(ax,"figures/L2/bottle/log/chla" + tmpT + ".png");
    clear tmp ax;

    % HPLC Monovinyl chlorophyll a (94-21)
    tmp = importdata("data\L2\mvchla_94-21_200.txt");
    ax = L2_helper(tmp,mld,dcm,thresh,testSel,"ks",[-60 60],[5 17],3);
    sgtitle("HPLC Monovinyl chlorophyll a (94-21): L2"+tmpT);
    exportgraphics(ax,"figures/L2/bottle/log/mvchla" + tmpT + ".png");
    clear tmp ax;
    
    % HPLC Divinyl chlorophyll a (94-21)
    tmp = importdata("data\L2\dvchla_94-21_200.txt");
    ax = L2_helper(tmp,mld,dcm,thresh,testSel,"ks",[-60 60],[5 17],3);
    sgtitle("HPLC Divinyl chlorophyll a (94-21): L2"+tmpT);
    exportgraphics(ax,"figures/L2/bottle/log/dvchla" + tmpT + ".png");
    clear tmp ax;
    
    % HPLC Chlorophyll b (88-21)
    tmp = importdata("data\L2\chlb_88-21_200.txt");
    ax = L2_helper(tmp,mld,dcm,thresh,testSel,"ks",[-60 60],[5 17],3);
    sgtitle("HPLC Chlorophyll b (88-21): L2"+tmpT);
    exportgraphics(ax,"figures/L2/bottle/log/chlb" + tmpT + ".png");
    clear tmp ax;
    
    % HPLC Chlorophyll C1 + C2 + C3 (88-21)
    tmp = importdata("data\L2\chl123_88-21_200.txt");
    ax = L2_helper(tmp,mld,dcm,thresh,testSel,"ks",[-60 60],[5 17],3);
    sgtitle("Chlorophyll C1 + C2 + C3 (88-21): L2"+tmpT);
    exportgraphics(ax,"figures/L2/bottle/log/chl123" + tmpT + ".png");
    clear tmp ax;
    
    % HPLC alpha-Carotene (94-21)
    tmp = importdata("data\L2\acar_94-21_200.txt");
    ax = L2_helper(tmp,mld,dcm,thresh,testSel,"ks",[-60 60],[5 17],3);
    sgtitle("HPLC alpha-Carotene (94-21): L2"+tmpT);
    exportgraphics(ax,"figures/L2/bottle/log/acar" + tmpT + ".png");
    clear tmp ax;
    
    % But19 (88-21)
    tmp = importdata("data\L2\but19_88-21_200.txt");
    ax = L2_helper(tmp,mld,dcm,thresh,testSel,"ks",[-60 60],[5 17],3);
    sgtitle("But19 (88-21): L2" + tmpT);
    exportgraphics(ax,"figures/L2/bottle/log/but19" + tmpT + ".png");
    clear tmp ax;
    
    % HPLC 19 Hexanoyloxyfucoxanthin (88-21)
    tmp = importdata("data\L2\hex19_88-21_200.txt");
    ax = L2_helper(tmp,mld,dcm,thresh,testSel,"ks",[-60 60],[5 17],3);
    sgtitle("HPLC 19' Hexanoyloxyfucoxanthin (88-21): L2" + tmpT);
    exportgraphics(ax,"figures/L2/bottle/log/hex19" + tmpT + ".png");
    clear tmp ax;
    
    % HPLC Zeaxanthin (88-21)
    tmp = importdata("data\L2\zeaxan_88-21_200.txt");
    ax = L2_helper(tmp,mld,dcm,thresh,testSel,"ks",[-60 60],[5 17],3);
    sgtitle("HPLC Zeaxanthin (88-21): L2" + tmpT);
    exportgraphics(ax,"figures/L2/bottle/log/zeaxan" + tmpT + ".png");
    clear tmp ax;
    
    % Particulate Carbon (89-21)
    tmp = importdata("data\L2\parc_89-21_200.txt");
    ax = L2_helper(tmp,mld,dcm,thresh,testSel,"ks",[-60 60],[5 17],3);
    sgtitle("Particulate Carbon 89-21: L2" + tmpT);
    exportgraphics(ax,"figures/L2/bottle/log/pc" + tmpT + ".png");
    clear tmp ax;
    
    % Particulate Nitrogen (89-21)
    tmp = importdata("data\L2\parn_89-21_200.txt");
    ax = L2_helper(tmp,mld,dcm,thresh,testSel,"ks",[-60 60],[5 17],3);
    sgtitle("Particulate Nitrogen 89-21: L2" + tmpT);
    exportgraphics(ax,"figures/L2/bottle/log/pn" + tmpT + ".png");
    clear tmp ax;
    
    % Nitrate + Nitrite (88-21)
    tmp = importdata("data/L2/nit2_88-21_200.txt");
    ax = L2_helper(tmp,mld,dcm,thresh,testSel,"ks",[-60 60],[5 17],3);
    sgtitle("Nitrate + Nitrite 88-21: L2" + tmpT);
    exportgraphics(ax,"figures/L2/bottle/log/nit2" + tmpT + ".png");
    clear tmp ax;
    
    % Phosphate (88-21)
    tmp = importdata("data/L2/pho_88-21_200.txt");
    ax = L2_helper(tmp,mld,dcm,thresh,testSel,"ks",[-60 60],[5 17],3);
    sgtitle("Phosphate 88-21: L2" + tmpT);
    exportgraphics(ax,"figures/L2/bottle/log/phos" + tmpT + ".png");
    clear tmp ax;
    
    % Oxygen (88-21)
    tmp = importdata("data/L2/oxy_88-21_200.txt");
    ax = L2_helper(tmp,mld,dcm,thresh,testSel,"ks",[-60 60],[6 18],3);
    sgtitle("Dissolved Oxygen 88-21: L2" + tmpT);
    exportgraphics(ax,"figures/L2/bottle/log/boxy" + tmpT + ".png");
    clear tmp ax;
    
    % PProd Light-12 (89-22)
    tmp = importdata("data\L2\l12_89-22_200.txt");
    ax = L2_helper(tmp,mld,dcm,10,testSel,"ks",[-60 60],[5 17],3);
    sgtitle("PProd Light-12 (89-22): L2" + tmpT);
    exportgraphics(ax,"figures/L2/bottle/log/l12" + tmpT + ".png");
    clear tmp ax;

    %%% AUTUMN
    tmpT = "-04";
    
    % chla
    tmp = importdata("data/L2/hplcChla_88-21_200.txt");
    ax = L2_helper(tmp,mld,dcm,15,testSel,"ks",[-50 50],[1 11],4);
    sgtitle("[Chl a] 88-21: L2"+tmpT);
    exportgraphics(ax,"figures/L2/bottle/log/chla" + tmpT + ".png");
    clear tmp ax;

    % HPLC Monovinyl chlorophyll a (94-21)
    tmp = importdata("data\L2\mvchla_94-21_200.txt");
    ax = L2_helper(tmp,mld,dcm,thresh,testSel,"ks",[-50 50],[1 11],4);
    sgtitle("HPLC Monovinyl chlorophyll a (94-21): L2"+tmpT);
    exportgraphics(ax,"figures/L2/bottle/log/mvchla" + tmpT + ".png");
    clear tmp ax;
    
    % HPLC Divinyl chlorophyll a (94-21)
    tmp = importdata("data\L2\dvchla_94-21_200.txt");
    ax = L2_helper(tmp,mld,dcm,thresh,testSel,"ks",[-50 50],[1 11],4);
    sgtitle("HPLC Divinyl chlorophyll a (94-21): L2"+tmpT);
    exportgraphics(ax,"figures/L2/bottle/log/dvchla" + tmpT + ".png");
    clear tmp ax;
    
    % HPLC Chlorophyll b (88-21)
    tmp = importdata("data\L2\chlb_88-21_200.txt");
    ax = L2_helper(tmp,mld,dcm,thresh,testSel,"ks",[-50 50],[1 11],4);
    sgtitle("HPLC Chlorophyll b (88-21): L2"+tmpT);
    exportgraphics(ax,"figures/L2/bottle/log/chlb" + tmpT + ".png");
    clear tmp ax;
    
    % HPLC Chlorophyll C1 + C2 + C3 (88-21)
    tmp = importdata("data\L2\chl123_88-21_200.txt");
    ax = L2_helper(tmp,mld,dcm,thresh,testSel,"ks",[-50 50],[1 11],4);
    sgtitle("Chlorophyll C1 + C2 + C3 (88-21): L2"+tmpT);
    exportgraphics(ax,"figures/L2/bottle/log/chl123" + tmpT + ".png");
    clear tmp ax;
    
    % HPLC alpha-Carotene (94-21)
    tmp = importdata("data\L2\acar_94-21_200.txt");
    ax = L2_helper(tmp,mld,dcm,thresh,testSel,"ks",[-50 50],[1 11],4);
    sgtitle("HPLC alpha-Carotene (94-21): L2"+tmpT);
    exportgraphics(ax,"figures/L2/bottle/log/acar" + tmpT + ".png");
    clear tmp ax;
    
    % But19 (88-21)
    tmp = importdata("data\L2\but19_88-21_200.txt");
    ax = L2_helper(tmp,mld,dcm,thresh,testSel,"ks",[-50 50],[1 11],4);
    sgtitle("But19 (88-21): L2" + tmpT);
    exportgraphics(ax,"figures/L2/bottle/log/but19" + tmpT + ".png");
    clear tmp ax;
    
    % HPLC 19 Hexanoyloxyfucoxanthin (88-21)
    tmp = importdata("data\L2\hex19_88-21_200.txt");
    ax = L2_helper(tmp,mld,dcm,thresh,testSel,"ks",[-50 50],[1 11],4);
    sgtitle("HPLC 19' Hexanoyloxyfucoxanthin (88-21): L2" + tmpT);
    exportgraphics(ax,"figures/L2/bottle/log/hex19" + tmpT + ".png");
    clear tmp ax;
    
    % HPLC Zeaxanthin (88-21)
    tmp = importdata("data\L2\zeaxan_88-21_200.txt");
    ax = L2_helper(tmp,mld,dcm,thresh,testSel,"ks",[-50 50],[1 11],4);
    sgtitle("HPLC Zeaxanthin (88-21): L2" + tmpT);
    exportgraphics(ax,"figures/L2/bottle/log/zeaxan" + tmpT + ".png");
    clear tmp ax;
    
    % Particulate Carbon (89-21)
    tmp = importdata("data\L2\parc_89-21_200.txt");
    ax = L2_helper(tmp,mld,dcm,thresh,testSel,"ks",[-50 50],[1 11],4);
    sgtitle("Particulate Carbon 89-21: L2" + tmpT);
    exportgraphics(ax,"figures/L2/bottle/log/pc" + tmpT + ".png");
    clear tmp ax;
    
    % Particulate Nitrogen (89-21)
    tmp = importdata("data\L2\parn_89-21_200.txt");
    ax = L2_helper(tmp,mld,dcm,thresh,testSel,"ks",[-50 50],[1 11],4);
    sgtitle("Particulate Nitrogen 89-21: L2" + tmpT);
    exportgraphics(ax,"figures/L2/bottle/log/pn" + tmpT + ".png");
    clear tmp ax;
    
    % Nitrate + Nitrite (88-21)
    tmp = importdata("data/L2/nit2_88-21_200.txt");
    ax = L2_helper(tmp,mld,dcm,thresh,testSel,"ks",[-60 60],[2 14],4);
    sgtitle("Nitrate + Nitrite 88-21: L2" + tmpT);
    exportgraphics(ax,"figures/L2/bottle/log/nit2" + tmpT + ".png");
    clear tmp ax;
    
    % Phosphate (88-21)
    tmp = importdata("data/L2/pho_88-21_200.txt");
    ax = L2_helper(tmp,mld,dcm,thresh,testSel,"ks",[-60 60],[2 14],4);
    sgtitle("Phosphate 88-21: L2" + tmpT);
    exportgraphics(ax,"figures/L2/bottle/log/phos" + tmpT + ".png");
    clear tmp ax;
    
    % Oxygen (88-21)
    tmp = importdata("data/L2/oxy_88-21_200.txt");
    ax = L2_helper(tmp,mld,dcm,thresh,testSel,"ks",[-50 50],[1 11],4);
    sgtitle("Dissolved Oxygen 88-21: L2" + tmpT);
    exportgraphics(ax,"figures/L2/bottle/log/boxy" + tmpT + ".png");
    clear tmp ax;
    
    % PProd Light-12 (89-22)
    tmp = importdata("data\L2\l12_89-22_200.txt");
    ax = L2_helper(tmp,mld,dcm,10,testSel,"ks",[-50 50],[1 11],4);
    sgtitle("PProd Light-12 (89-22): L2" + tmpT);
    exportgraphics(ax,"figures/L2/bottle/log/l12" + tmpT + ".png");
    clear tmp ax;
  
end

%% Seasonal Analysis: A-D

if seasonalAnalysisAd == true

    set(groot, "defaultFigureUnits", "centimeters", "defaultFigurePosition", [3 3 10 15]);

    thresh = 30;
    %%% WINTER
    tmpT = "-ad-01";

    % Chlorophyll a (88-21)
    tmpx = "Winter";
    tmp = importdata("data/L2/hplcChla_88-21_200.txt");
    ax = L2_helper(tmp,mld,dcm,thresh,testSel,"ad",[-50 50],[4 14],1);
    sgtitle("L2 "+tmpx,"Interpreter","latex");
    exportgraphics(ax,"figures/L2/bottle/log/chla" + tmpT + ".png");
    clear tmp ax;

    % HPLC Monovinyl chlorophyll a (94-21)
    tmp = importdata("data\L2\mvchla_94-21_200.txt");
    ax = L2_helper(tmp,mld,dcm,thresh,testSel,"ad",[-60 60],[5 17],1);
    sgtitle("HPLC Monovinyl chlorophyll a (94-21): L2"+tmpT);
    exportgraphics(ax,"figures/L2/bottle/log/mvchla" + tmpT + ".png");
    clear tmp ax;
    
    % HPLC Divinyl chlorophyll a (94-21)
    tmp = importdata("data\L2\dvchla_94-21_200.txt");
    ax = L2_helper(tmp,mld,dcm,thresh,testSel,"ad",[-60 60],[5 17],1);
    sgtitle("HPLC Divinyl chlorophyll a (94-21): L2"+tmpT);
    exportgraphics(ax,"figures/L2/bottle/log/dvchla" + tmpT + ".png");
    clear tmp ax;
    
    % HPLC Chlorophyll b (88-21)
    tmp = importdata("data\L2\chlb_88-21_200.txt");
    ax = L2_helper(tmp,mld,dcm,thresh,testSel,"ad",[-60 60],[5 17],1);
    sgtitle("HPLC Chlorophyll b (88-21): L2"+tmpT);
    exportgraphics(ax,"figures/L2/bottle/log/chlb" + tmpT + ".png");
    clear tmp ax;
    
    % HPLC Chlorophyll C1 + C2 + C3 (88-21)
    tmp = importdata("data\L2\chl123_88-21_200.txt");
    ax = L2_helper(tmp,mld,dcm,thresh,testSel,"ad",[-60 60],[5 17],1);
    sgtitle("Chlorophyll C1 + C2 + C3 (88-21): L2"+tmpT);
    exportgraphics(ax,"figures/L2/bottle/log/chl123" + tmpT + ".png");
    clear tmp ax;
    
    % HPLC alpha-Carotene (94-21)
    tmp = importdata("data\L2\acar_94-21_200.txt");
    ax = L2_helper(tmp,mld,dcm,thresh,testSel,"ad",[-60 60],[5 17],1);
    sgtitle("HPLC alpha-Carotene (94-21): L2"+tmpT);
    exportgraphics(ax,"figures/L2/bottle/log/acar" + tmpT + ".png");
    clear tmp ax;
    
    % But19 (88-21)
    tmp = importdata("data\L2\but19_88-21_200.txt");
    ax = L2_helper(tmp,mld,dcm,thresh,testSel,"ad",[-60 60],[5 17],1);
    sgtitle("But19 (88-21): L2" + tmpT);
    exportgraphics(ax,"figures/L2/bottle/log/but19" + tmpT + ".png");
    clear tmp ax;
    
    % HPLC 19 Hexanoyloxyfucoxanthin (88-21)
    tmp = importdata("data\L2\hex19_88-21_200.txt");
    ax = L2_helper(tmp,mld,dcm,thresh,testSel,"ad",[-60 60],[5 17],1);
    sgtitle("HPLC 19' Hexanoyloxyfucoxanthin (88-21): L2" + tmpT);
    exportgraphics(ax,"figures/L2/bottle/log/hex19" + tmpT + ".png");
    clear tmp ax;
    
    % HPLC Zeaxanthin (88-21)
    tmp = importdata("data\L2\zeaxan_88-21_200.txt");
    ax = L2_helper(tmp,mld,dcm,thresh,testSel,"ad",[-60 60],[5 17],1);
    sgtitle("HPLC Zeaxanthin (88-21): L2" + tmpT);
    exportgraphics(ax,"figures/L2/bottle/log/zeaxan" + tmpT + ".png");
    clear tmp ax;
    
    % Particulate Carbon (89-21)
    tmp = importdata("data\L2\parc_89-21_200.txt");
    ax = L2_helper(tmp,mld,dcm,thresh,testSel,"ad",[-60 60],[1 13],1);
    sgtitle("L2 Winter");
    exportgraphics(ax,"figures/L2/bottle/log/pc" + tmpT + ".png");
    clear tmp ax;
    
    % Particulate Nitrogen (89-21)
    tmp = importdata("data\L2\parn_89-21_200.txt");
    ax = L2_helper(tmp,mld,dcm,thresh,testSel,"ad",[-60 60],[2 14],1);
    sgtitle("Particulate Nitrogen 89-21: L2" + tmpT);
    exportgraphics(ax,"figures/L2/bottle/log/pn" + tmpT + ".png");
    clear tmp ax;
    
    % Nitrate + Nitrite (88-21)
    tmp = importdata("data/L2/nit2_88-21_200.txt");
    ax = L2_helper(tmp,mld,dcm,thresh,testSel,"ad",[-60 60],[6 18],1);
    sgtitle("Nitrate + Nitrite 88-21: L2" + tmpT);
    exportgraphics(ax,"figures/L2/bottle/log/nit2" + tmpT + ".png");
    clear tmp ax;
    
    % Phosphate (88-21)
    tmp = importdata("data/L2/pho_88-21_200.txt");
    ax = L2_helper(tmp,mld,dcm,thresh,testSel,"ad",[-60 60],[6 18],1);
    sgtitle("Phosphate 88-21: L2" + tmpT);
    exportgraphics(ax,"figures/L2/bottle/log/phos" + tmpT + ".png");
    clear tmp ax;
    
    % Oxygen (88-21)
    tmp = importdata("data/L2/oxy_88-21_200.txt");
    ax = L2_helper(tmp,mld,dcm,thresh,testSel,"ad",[-60 60],[5 17],1);
    sgtitle("Dissolved Oxygen 88-21: L2" + tmpT);
    exportgraphics(ax,"figures/L2/bottle/log/boxy" + tmpT + ".png");
    clear tmp ax;
    
    % PProd Light-12 (89-22)
    tmp = importdata("data\L2\l12_89-22_200.txt");
    ax = L2_helper(tmp,mld,dcm,thresh,testSel,"ad",[-60 60],[3 15],1);
    sgtitle("PProd Light-12 (89-22): L2" + tmpT);
    exportgraphics(ax,"figures/L2/bottle/log/l12" + tmpT + ".png");
    clear tmp ax;


    %%% SPRING
    tmpT = "-ad-02";

    % chla
    tmpx = "Spring";
    tmp = importdata("data/L2/hplcChla_88-21_200.txt");
    ax = L2_helper(tmp,mld,dcm,thresh,testSel,"ad",[-50 50],[7 17],2);
    sgtitle("L2 "+tmpx,"Interpreter","latex");
    exportgraphics(ax,"figures/L2/bottle/log/chla" + tmpT + ".png");
    clear tmp ax;

    % HPLC Monovinyl chlorophyll a (94-21)
    tmp = importdata("data\L2\mvchla_94-21_200.txt");
    ax = L2_helper(tmp,mld,dcm,thresh,testSel,"ad",[-60 60],[7 19],2);
    sgtitle("HPLC Monovinyl chlorophyll a (94-21): L2"+tmpT);
    exportgraphics(ax,"figures/L2/bottle/log/mvchla" + tmpT + ".png");
    clear tmp ax;
    
    % HPLC Divinyl chlorophyll a (94-21)
    tmp = importdata("data\L2\dvchla_94-21_200.txt");
    ax = L2_helper(tmp,mld,dcm,thresh,testSel,"ad",[-60 60],[7 19],2);
    sgtitle("HPLC Divinyl chlorophyll a (94-21): L2"+tmpT);
    exportgraphics(ax,"figures/L2/bottle/log/dvchla" + tmpT + ".png");
    clear tmp ax;
    
    % HPLC Chlorophyll b (88-21)
    tmp = importdata("data\L2\chlb_88-21_200.txt");
    ax = L2_helper(tmp,mld,dcm,thresh,testSel,"ad",[-60 60],[7 19],2);
    sgtitle("HPLC Chlorophyll b (88-21): L2"+tmpT);
    exportgraphics(ax,"figures/L2/bottle/log/chlb" + tmpT + ".png");
    clear tmp ax;
    
    % HPLC Chlorophyll C1 + C2 + C3 (88-21)
    tmp = importdata("data\L2\chl123_88-21_200.txt");
    ax = L2_helper(tmp,mld,dcm,thresh,testSel,"ad",[-60 60],[7 19],2);
    sgtitle("Chlorophyll C1 + C2 + C3 (88-21): L2"+tmpT);
    exportgraphics(ax,"figures/L2/bottle/log/chl123" + tmpT + ".png");
    clear tmp ax;
    
    % HPLC alpha-Carotene (94-21)
    tmp = importdata("data\L2\acar_94-21_200.txt");
    ax = L2_helper(tmp,mld,dcm,thresh,testSel,"ad",[-60 60],[7 19],2);
    sgtitle("HPLC alpha-Carotene (94-21): L2"+tmpT);
    exportgraphics(ax,"figures/L2/bottle/log/acar" + tmpT + ".png");
    clear tmp ax;
    
    % But19 (88-21)
    tmp = importdata("data\L2\but19_88-21_200.txt");
    ax = L2_helper(tmp,mld,dcm,thresh,testSel,"ad",[-60 60],[7 19],2);
    sgtitle("But19 (88-21): L2" + tmpT);
    exportgraphics(ax,"figures/L2/bottle/log/but19" + tmpT + ".png");
    clear tmp ax;
    
    % HPLC 19 Hexanoyloxyfucoxanthin (88-21)
    tmp = importdata("data\L2\hex19_88-21_200.txt");
    ax = L2_helper(tmp,mld,dcm,thresh,testSel,"ad",[-60 60],[7 19],2);
    sgtitle("HPLC 19' Hexanoyloxyfucoxanthin (88-21): L2" + tmpT);
    exportgraphics(ax,"figures/L2/bottle/log/hex19" + tmpT + ".png");
    clear tmp ax;
    
    % HPLC Zeaxanthin (88-21)
    tmp = importdata("data\L2\zeaxan_88-21_200.txt");
    ax = L2_helper(tmp,mld,dcm,thresh,testSel,"ad",[-60 60],[7 19],2);
    sgtitle("HPLC Zeaxanthin (88-21): L2" + tmpT);
    exportgraphics(ax,"figures/L2/bottle/log/zeaxan" + tmpT + ".png");
    clear tmp ax;
    
    % Particulate Carbon (89-21)
    tmp = importdata("data\L2\parc_89-21_200.txt");
    ax = L2_helper(tmp,mld,dcm,thresh,testSel,"ad",[-60 60],[8 20],2);
    sgtitle("L2 Spring");
    exportgraphics(ax,"figures/L2/bottle/log/pc" + tmpT + ".png");
    clear tmp ax;
    
    % Particulate Nitrogen (89-21)
    tmp = importdata("data\L2\parn_89-21_200.txt");
    ax = L2_helper(tmp,mld,dcm,thresh,testSel,"ad",[-60 60],[8 20],2);
    sgtitle("Particulate Nitrogen 89-21: L2" + tmpT);
    exportgraphics(ax,"figures/L2/bottle/log/pn" + tmpT + ".png");
    clear tmp ax;
    
    % Nitrate + Nitrite (88-21)
    tmp = importdata("data/L2/nit2_88-21_200.txt");
    ax = L2_helper(tmp,mld,dcm,thresh,testSel,"ad",[-60 60],[6 18],2);
    sgtitle("Nitrate + Nitrite 88-21: L2" + tmpT);
    exportgraphics(ax,"figures/L2/bottle/log/nit2" + tmpT + ".png");
    clear tmp ax;
    
    % Phosphate (88-21)
    tmp = importdata("data/L2/pho_88-21_200.txt");
    ax = L2_helper(tmp,mld,dcm,thresh,testSel,"ad",[-60 60],[6 18],2);
    sgtitle("Phosphate 88-21: L2" + tmpT);
    exportgraphics(ax,"figures/L2/bottle/log/phos" + tmpT + ".png");
    clear tmp ax;
    
    % Oxygen (88-21)
    tmp = importdata("data/L2/oxy_88-21_200.txt");
    ax = L2_helper(tmp,mld,dcm,thresh,testSel,"ad",[-60 60],[6 18],2);
    sgtitle("Dissolved Oxygen 88-21: L2" + tmpT);
    exportgraphics(ax,"figures/L2/bottle/log/boxy" + tmpT + ".png");
    clear tmp ax;
    
    % PProd Light-12 (89-22)
    tmp = importdata("data\L2\l12_89-22_200.txt");
    ax = L2_helper(tmp,mld,dcm,10,testSel,"ad",[-60 60],[7 19],2);
    sgtitle("PProd Light-12 (89-22): L2" + tmpT);
    exportgraphics(ax,"figures/L2/bottle/log/l12" + tmpT + ".png");
    clear tmp ax;


    %%% SUMMER
    tmpT = "-ad-03";

    % chla
    tmpx = " Summer";
    tmp = importdata("data/L2/hplcChla_88-21_200.txt");
    ax = L2_helper(tmp,mld,dcm,thresh,testSel,"ad",[-50 50],[6 16],3);
    sgtitle("L2"+tmpx,"Interpreter","latex");
    exportgraphics(ax,"figures/L2/bottle/log/chla" + tmpT + ".png");
    clear tmp ax;

    % HPLC Monovinyl chlorophyll a (94-21)
    tmp = importdata("data\L2\mvchla_94-21_200.txt");
    ax = L2_helper(tmp,mld,dcm,thresh,testSel,"ad",[-60 60],[5 17],3);
    sgtitle("HPLC Monovinyl chlorophyll a (94-21): L2"+tmpT);
    exportgraphics(ax,"figures/L2/bottle/log/mvchla" + tmpT + ".png");
    clear tmp ax;
    
    % HPLC Divinyl chlorophyll a (94-21)
    tmp = importdata("data\L2\dvchla_94-21_200.txt");
    ax = L2_helper(tmp,mld,dcm,thresh,testSel,"ad",[-60 60],[5 17],3);
    sgtitle("HPLC Divinyl chlorophyll a (94-21): L2"+tmpT);
    exportgraphics(ax,"figures/L2/bottle/log/dvchla" + tmpT + ".png");
    clear tmp ax;
    
    % HPLC Chlorophyll b (88-21)
    tmp = importdata("data\L2\chlb_88-21_200.txt");
    ax = L2_helper(tmp,mld,dcm,thresh,testSel,"ad",[-60 60],[5 17],3);
    sgtitle("HPLC Chlorophyll b (88-21): L2"+tmpT);
    exportgraphics(ax,"figures/L2/bottle/log/chlb" + tmpT + ".png");
    clear tmp ax;
    
    % HPLC Chlorophyll C1 + C2 + C3 (88-21)
    tmp = importdata("data\L2\chl123_88-21_200.txt");
    ax = L2_helper(tmp,mld,dcm,thresh,testSel,"ad",[-60 60],[5 17],3);
    sgtitle("Chlorophyll C1 + C2 + C3 (88-21): L2"+tmpT);
    exportgraphics(ax,"figures/L2/bottle/log/chl123" + tmpT + ".png");
    clear tmp ax;
    
    % HPLC alpha-Carotene (94-21)
    tmp = importdata("data\L2\acar_94-21_200.txt");
    ax = L2_helper(tmp,mld,dcm,thresh,testSel,"ad",[-60 60],[5 17],3);
    sgtitle("HPLC alpha-Carotene (94-21): L2"+tmpT);
    exportgraphics(ax,"figures/L2/bottle/log/acar" + tmpT + ".png");
    clear tmp ax;
    
    % But19 (88-21)
    tmp = importdata("data\L2\but19_88-21_200.txt");
    ax = L2_helper(tmp,mld,dcm,thresh,testSel,"ad",[-60 60],[5 17],3);
    sgtitle("But19 (88-21): L2" + tmpT);
    exportgraphics(ax,"figures/L2/bottle/log/but19" + tmpT + ".png");
    clear tmp ax;
    
    % HPLC 19 Hexanoyloxyfucoxanthin (88-21)
    tmp = importdata("data\L2\hex19_88-21_200.txt");
    ax = L2_helper(tmp,mld,dcm,thresh,testSel,"ad",[-60 60],[5 17],3);
    sgtitle("HPLC 19' Hexanoyloxyfucoxanthin (88-21): L2" + tmpT);
    exportgraphics(ax,"figures/L2/bottle/log/hex19" + tmpT + ".png");
    clear tmp ax;
    
    % HPLC Zeaxanthin (88-21)
    tmp = importdata("data\L2\zeaxan_88-21_200.txt");
    ax = L2_helper(tmp,mld,dcm,thresh,testSel,"ad",[-60 60],[5 17],3);
    sgtitle("HPLC Zeaxanthin (88-21): L2" + tmpT);
    exportgraphics(ax,"figures/L2/bottle/log/zeaxan" + tmpT + ".png");
    clear tmp ax;
    
    % Particulate Carbon (89-21)
    tmp = importdata("data\L2\parc_89-21_200.txt");
    ax = L2_helper(tmp,mld,dcm,thresh,testSel,"ad",[-60 60],[5 17],3);
    sgtitle("L2 Summer");
    exportgraphics(ax,"figures/L2/bottle/log/pc" + tmpT + ".png");
    clear tmp ax;
    
    % Particulate Nitrogen (89-21)
    tmp = importdata("data\L2\parn_89-21_200.txt");
    ax = L2_helper(tmp,mld,dcm,thresh,testSel,"ad",[-60 60],[5 17],3);
    sgtitle("Particulate Nitrogen 89-21: L2" + tmpT);
    exportgraphics(ax,"figures/L2/bottle/log/pn" + tmpT + ".png");
    clear tmp ax;
    
    % Nitrate + Nitrite (88-21)
    tmp = importdata("data/L2/nit2_88-21_200.txt");
    ax = L2_helper(tmp,mld,dcm,thresh,testSel,"ad",[-60 60],[5 17],3);
    sgtitle("Nitrate + Nitrite 88-21: L2" + tmpT);
    exportgraphics(ax,"figures/L2/bottle/log/nit2" + tmpT + ".png");
    clear tmp ax;
    
    % Phosphate (88-21)
    tmp = importdata("data/L2/pho_88-21_200.txt");
    ax = L2_helper(tmp,mld,dcm,thresh,testSel,"ad",[-60 60],[5 17],3);
    sgtitle("Phosphate 88-21: L2" + tmpT);
    exportgraphics(ax,"figures/L2/bottle/log/phos" + tmpT + ".png");
    clear tmp ax;
    
    % Oxygen (88-21)
    tmp = importdata("data/L2/oxy_88-21_200.txt");
    ax = L2_helper(tmp,mld,dcm,thresh,testSel,"ad",[-60 60],[6 18],3);
    sgtitle("Dissolved Oxygen 88-21: L2" + tmpT);
    exportgraphics(ax,"figures/L2/bottle/log/boxy" + tmpT + ".png");
    clear tmp ax;
    
    % PProd Light-12 (89-22)
    tmp = importdata("data\L2\l12_89-22_200.txt");
    ax = L2_helper(tmp,mld,dcm,10,testSel,"ad",[-60 60],[5 17],3);
    sgtitle("PProd Light-12 (89-22): L2" + tmpT);
    exportgraphics(ax,"figures/L2/bottle/log/l12" + tmpT + ".png");
    clear tmp ax;


    %%% AUTUMN
    tmpT = "-ad-04";
    
    % chla
    tmpx = " Autumn";
    tmp = importdata("data/L2/hplcChla_88-21_200.txt");
    ax = L2_helper(tmp,mld,dcm,thresh,testSel,"ad",[-50 50],[3 13],4);
    sgtitle("L2"+tmpx,"Interpreter","latex");
    exportgraphics(ax,"figures/L2/bottle/log/chla" + tmpT + ".png");
    clear tmp ax;

    % HPLC Monovinyl chlorophyll a (94-21)
    tmp = importdata("data\L2\mvchla_94-21_200.txt");
    ax = L2_helper(tmp,mld,dcm,thresh,testSel,"ad",[-50 50],[1 11],4);
    sgtitle("HPLC Monovinyl chlorophyll a (94-21): L2"+tmpT);
    exportgraphics(ax,"figures/L2/bottle/log/mvchla" + tmpT + ".png");
    clear tmp ax;
    
    % HPLC Divinyl chlorophyll a (94-21)
    tmp = importdata("data\L2\dvchla_94-21_200.txt");
    ax = L2_helper(tmp,mld,dcm,thresh,testSel,"ad",[-50 50],[1 11],4);
    sgtitle("HPLC Divinyl chlorophyll a (94-21): L2"+tmpT);
    exportgraphics(ax,"figures/L2/bottle/log/dvchla" + tmpT + ".png");
    clear tmp ax;
    
    % HPLC Chlorophyll b (88-21)
    tmp = importdata("data\L2\chlb_88-21_200.txt");
    ax = L2_helper(tmp,mld,dcm,thresh,testSel,"ad",[-50 50],[1 11],4);
    sgtitle("HPLC Chlorophyll b (88-21): L2"+tmpT);
    exportgraphics(ax,"figures/L2/bottle/log/chlb" + tmpT + ".png");
    clear tmp ax;
    
    % HPLC Chlorophyll C1 + C2 + C3 (88-21)
    tmp = importdata("data\L2\chl123_88-21_200.txt");
    ax = L2_helper(tmp,mld,dcm,thresh,testSel,"ad",[-50 50],[1 11],4);
    sgtitle("Chlorophyll C1 + C2 + C3 (88-21): L2"+tmpT);
    exportgraphics(ax,"figures/L2/bottle/log/chl123" + tmpT + ".png");
    clear tmp ax;
    
    % HPLC alpha-Carotene (94-21)
    tmp = importdata("data\L2\acar_94-21_200.txt");
    ax = L2_helper(tmp,mld,dcm,thresh,testSel,"ad",[-50 50],[1 11],4);
    sgtitle("HPLC alpha-Carotene (94-21): L2"+tmpT);
    exportgraphics(ax,"figures/L2/bottle/log/acar" + tmpT + ".png");
    clear tmp ax;
    
    % But19 (88-21)
    tmp = importdata("data\L2\but19_88-21_200.txt");
    ax = L2_helper(tmp,mld,dcm,thresh,testSel,"ad",[-50 50],[1 11],4);
    sgtitle("But19 (88-21): L2" + tmpT);
    exportgraphics(ax,"figures/L2/bottle/log/but19" + tmpT + ".png");
    clear tmp ax;
    
    % HPLC 19 Hexanoyloxyfucoxanthin (88-21)
    tmp = importdata("data\L2\hex19_88-21_200.txt");
    ax = L2_helper(tmp,mld,dcm,thresh,testSel,"ad",[-50 50],[1 11],4);
    sgtitle("HPLC 19' Hexanoyloxyfucoxanthin (88-21): L2" + tmpT);
    exportgraphics(ax,"figures/L2/bottle/log/hex19" + tmpT + ".png");
    clear tmp ax;
    
    % HPLC Zeaxanthin (88-21)
    tmp = importdata("data\L2\zeaxan_88-21_200.txt");
    ax = L2_helper(tmp,mld,dcm,thresh,testSel,"ad",[-50 50],[1 11],4);
    sgtitle("HPLC Zeaxanthin (88-21): L2" + tmpT);
    exportgraphics(ax,"figures/L2/bottle/log/zeaxan" + tmpT + ".png");
    clear tmp ax;
    
    % Particulate Carbon (89-21)
    tmp = importdata("data\L2\parc_89-21_200.txt");
    ax = L2_helper(tmp,mld,dcm,thresh,testSel,"ad",[-50 50],[1 11],4);
    sgtitle("L2 Autumn");
    exportgraphics(ax,"figures/L2/bottle/log/pc" + tmpT + ".png");
    clear tmp ax;
    
    % Particulate Nitrogen (89-21)
    tmp = importdata("data\L2\parn_89-21_200.txt");
    ax = L2_helper(tmp,mld,dcm,thresh,testSel,"ad",[-50 50],[1 11],4);
    sgtitle("Particulate Nitrogen 89-21: L2" + tmpT);
    exportgraphics(ax,"figures/L2/bottle/log/pn" + tmpT + ".png");
    clear tmp ax;
    
    % Nitrate + Nitrite (88-21)
    tmp = importdata("data/L2/nit2_88-21_200.txt");
    ax = L2_helper(tmp,mld,dcm,thresh,testSel,"ad",[-60 60],[2 14],4);
    sgtitle("Nitrate + Nitrite 88-21: L2" + tmpT);
    exportgraphics(ax,"figures/L2/bottle/log/nit2" + tmpT + ".png");
    clear tmp ax;
    
    % Phosphate (88-21)
    tmp = importdata("data/L2/pho_88-21_200.txt");
    ax = L2_helper(tmp,mld,dcm,thresh,testSel,"ad",[-60 60],[2 14],4);
    sgtitle("Phosphate 88-21: L2" + tmpT);
    exportgraphics(ax,"figures/L2/bottle/log/phos" + tmpT + ".png");
    clear tmp ax;
    
    % Oxygen (88-21)
    tmp = importdata("data/L2/oxy_88-21_200.txt");
    ax = L2_helper(tmp,mld,dcm,thresh,testSel,"ad",[-50 50],[1 11],4);
    sgtitle("Dissolved Oxygen 88-21: L2" + tmpT);
    exportgraphics(ax,"figures/L2/bottle/log/boxy" + tmpT + ".png");
    clear tmp ax;
    
    % PProd Light-12 (89-22)
    tmp = importdata("data\L2\l12_89-22_200.txt");
    ax = L2_helper(tmp,mld,dcm,10,testSel,"ad",[-50 50],[1 11],4);
    sgtitle("PProd Light-12 (89-22): L2" + tmpT);
    exportgraphics(ax,"figures/L2/bottle/log/l12" + tmpT + ".png");
    clear tmp ax;
end

%% Principal Analysis
if principleAnalysis == true
    %  K-S
    tmpT = "";
    
    % Chlorophyll a (88-21)
    tmp = importdata("data/L2/hplcChla_88-21_200.txt");
    [ax,p,ks,obs,sk,ku,~,~,~,ad,pr,VuongRes] = L2_helper(tmp,mld,dcm,50,testSel,"ks",[-60 60],[7 19]);
    sgtitle("[Chl a] 88-21: L2");
    exportgraphics(ax,"figures/L2/bottle/log/chla" + tmpT + ".png");
    save("output\L2\chla.mat","p","ks","obs","sk","ku");
    clearvars -except mld dcm tmpT testSel;
    
    % HPLC Monovinyl chlorophyll a (94-21)
    tmp = importdata("data\L2\mvchla_94-21_200.txt");
    [ax,p,ks,obs,sk,ku] = L2_helper(tmp,mld,dcm,50,testSel,"ks",[-60 60],[7 19]);
    sgtitle("HPLC Monovinyl chlorophyll a (94-21): L2");
    exportgraphics(ax,"figures/L2/bottle/log/mvchla" + tmpT + ".png");
    save("output\L2\mvchla.mat","p","ks","obs","sk","ku");
    clearvars -except mld dcm tmpT testSel;
    
    % HPLC Divinyl chlorophyll a (94-21)
    tmp = importdata("data\L2\dvchla_94-21_200.txt");
    [ax,p,ks,obs,sk,ku] = L2_helper(tmp,mld,dcm,50,testSel,"ks",[-60 60],[7 19]);
    sgtitle("HPLC Divinyl chlorophyll a (94-21): L2");
    exportgraphics(ax,"figures/L2/bottle/log/dvchla" + tmpT + ".png");
    save("output\L2\dvchla.mat","p","ks","obs","sk","ku");
    clearvars -except mld dcm tmpT testSel;
    
    % HPLC Chlorophyll b (88-21)
    tmp = importdata("data\L2\chlb_88-21_200.txt");
    [ax,p,ks,obs,sk,ku] = L2_helper(tmp,mld,dcm,50,testSel,"ks",[-60 60],[7 19]);
    sgtitle("HPLC Chlorophyll b (88-21): L2");
    exportgraphics(ax,"figures/L2/bottle/log/chlb" + tmpT + ".png");
    save("output\L2\chlb.mat","p","ks","obs","sk","ku");
    clearvars -except mld dcm tmpT testSel;
    
    % HPLC Chlorophyll C1 + C2 + C3 (88-21)
    tmp = importdata("data\L2\chl123_88-21_200.txt");
    [ax,p,ks,obs,sk,ku] = L2_helper(tmp,mld,dcm,50,testSel,"ks",[-60 60],[7 19]);
    sgtitle("Chlorophyll C1 + C2 + C3 (88-21): L2");
    exportgraphics(ax,"figures/L2/bottle/log/chl123" + tmpT + ".png");
    save("output\L2\chl123.mat","p","ks","obs","sk","ku");
    clearvars -except mld dcm tmpT testSel;
    
    % HPLC alpha-Carotene (94-21)
    tmp = importdata("data\L2\acar_94-21_200.txt");
    [ax,p,ks,obs,sk,ku] = L2_helper(tmp,mld,dcm,50,testSel,"ks",[-60 60],[7 19]);
    sgtitle("HPLC alpha-Carotene (94-21): L2");
    exportgraphics(ax,"figures/L2/bottle/log/acar" + tmpT + ".png");
    save("output\L2\acar.mat","p","ks","obs","sk","ku");
    clearvars -except mld dcm tmpT testSel;
    
    % But19 (88-21)
    tmp = importdata("data\L2\but19_88-21_200.txt");
    [ax,p,ks,obs,sk,ku] = L2_helper(tmp,mld,dcm,50,testSel,"ks",[-60 60],[7 19]);
    sgtitle("But19 (88-21): L2");
    exportgraphics(ax,"figures/L2/bottle/log/but19" + tmpT + ".png");
    save("output\L2\but19.mat","p","ks","obs","sk","ku");
    clearvars -except mld dcm tmpT testSel;
    
    % HPLC 19 Hexanoyloxyfucoxanthin (88-21)
    tmp = importdata("data\L2\hex19_88-21_200.txt");
    [ax,p,ks,obs,sk,ku] = L2_helper(tmp,mld,dcm,50,testSel,"ks",[-60 60],[7 19]);
    sgtitle("HPLC 19' Hexanoyloxyfucoxanthin (88-21): L2");
    exportgraphics(ax,"figures/L2/bottle/log/hex19" + tmpT + ".png");
    save("output\L2\hex19.mat","p","ks","obs","sk","ku");
    clearvars -except mld dcm tmpT testSel;
    
    % HPLC Zeaxanthin (88-21)
    tmp = importdata("data\L2\zeaxan_88-21_200.txt");
    [ax,p,ks,obs,sk,ku] = L2_helper(tmp,mld,dcm,50,testSel,"ks",[-60 60],[7 19]);
    sgtitle("HPLC Zeaxanthin (88-21): L2");
    exportgraphics(ax,"figures/L2/bottle/log/zeaxan" + tmpT + ".png");
    save("output\L2\zeaxan.mat","p","ks","obs","sk","ku");
    clearvars -except mld dcm tmpT testSel;
    
    % Particulate Carbon (89-21)
    tmp = importdata("data\L2\parc_89-21_200.txt");
    [ax,p,ks,obs,sk,ku] = L2_helper(tmp,mld,dcm,50,testSel,"ks",[-80 80],[6 22]);
    sgtitle("Particulate Carbon 89-21: L2");
    exportgraphics(ax,"figures/L2/bottle/log/pc" + tmpT + ".png");
    save("output\L2\pc.mat","p","ks","obs","sk","ku");
    clearvars -except mld dcm tmpT testSel;
    
    % Particulate Nitrogen (89-21)
    tmp = importdata("data\L2\parn_89-21_200.txt");
    [ax,p,ks,obs,sk,ku] = L2_helper(tmp,mld,dcm,50,testSel,"ks",[-60 60],[8 20]);
    sgtitle("Particulate Nitrogen 89-21: L2");
    exportgraphics(ax,"figures/L2/bottle/log/pn" + tmpT + ".png");
    save("output\L2\pn.mat","p","ks","obs","sk","ku");
    clearvars -except mld dcm tmpT testSel;
    
    % Nitrate + Nitrite (88-21)
    tmp = importdata("data/L2/nit2_88-21_200.txt");
    [ax,p,ks,obs,sk,ku] = L2_helper(tmp,mld,dcm,50,testSel,"ks",[-60 60],[7 19]);
    sgtitle("Nitrate + Nitrite 88-21: L2");
    exportgraphics(ax,"figures/L2/bottle/log/nit2" + tmpT + ".png");
    save("output\L2\nit2.mat","p","ks","obs","sk","ku");
    clearvars -except mld dcm tmpT testSel;
    
    % Phosphate (88-21)
    tmp = importdata("data/L2/pho_88-21_200.txt");
    [ax,p,ks,obs,sk,ku] = L2_helper(tmp,mld,dcm,50,testSel,"ks",[-60 60],[7 19]);
    sgtitle("Phosphate 88-21: L2");
    exportgraphics(ax,"figures/L2/bottle/log/phos" + tmpT + ".png");
    save("output\L2\phos.mat","p","ks","obs","sk","ku");
    clearvars -except mld dcm tmpT testSel;
    
    % Oxygen (88-21)
    tmp = importdata("data/L2/oxy_88-21_200.txt");
    [ax,p,ks,obs,sk,ku,rV,pSubml,pV,ad,pr,vuongRes] = L2_helper(tmp,mld,dcm,50,testSel,"ks",[-60 60],[6 18]);
    sgtitle("Dissolved Oxygen 88-21: L2");
    exportgraphics(ax,"figures/L2/bottle/log/boxy" + tmpT + ".png");
    save("output\L2\boxy.mat","p","ks","obs","sk","ku");
    clearvars -except mld dcm tmpT testSel;
    
    % PProd Light-12 (89-22)
    tmp = importdata("data\L2\l12_89-22_200.txt");
    [ax,p,ks,obs,sk,ku] = L2_helper(tmp,mld,dcm,50,testSel,"ks",[-60 60],[6 18]);
    sgtitle("PProd Light-12 (89-22): L2");
    exportgraphics(ax,"figures/L2/bottle/log/l12" + tmpT + ".png");
    save("output\L2\l12.mat","p","ks","obs","sk","ku");
    clearvars -except mld dcm tmpT testSel;

%     % Macrozooplankton (94-22)
%     tmp = importdata("data/L0/macrozoo_94-22_200.txt");
%     ax = L2_helper(tmp,mld,dcm,50,testSel,"ks",[-60 60],[7 19]);
%     sgtitle("Macrozooplankton 94-21: L2");
%     exportgraphics(ax,"figures/L2/bottle/log/macrozoo" + tmpT + ".png");
%     %save("output\L2\chla.mat","p","ks","obs","sk","ku");
%     %clearvars -except mld dcm tmpT testSel;
    
    % A-D
    tmpT = "-ad";
    
    testSel = 2;
    
    % Chlorophyll a (88-21)
    tmpX = "";
    tmp = importdata("data/L2/hplcChla_88-21_200.txt");
    ax = L2_helper(tmp,mld,dcm,30,testSel,"ad",[-60 60],[7 19]);
    sgtitle("L2"+tmpX,"Interpreter","latex");
    exportgraphics(ax,"figures/L2/bottle/log/chla" + tmpT + ".png");
    % save("output\L2\chla.mat","p","ks","obs","sk","ku");
    clearvars -except mld dcm tmpT testSel;
    
    % HPLC Monovinyl chlorophyll a (94-21)
    tmp = importdata("data\L2\mvchla_94-21_200.txt");
    ax = L2_helper(tmp,mld,dcm,50,testSel,"ad",[-60 60],[7 19]);
    sgtitle("HPLC Monovinyl chlorophyll a (94-21): L2"+tmpT);
    exportgraphics(ax,"figures/L2/bottle/log/mvchla" + tmpT + ".png");
    % save("output\L2\mvchla.mat","p","ks","obs","sk","ku");
    clearvars -except mld dcm tmpT testSel;
    
    % HPLC Divinyl chlorophyll a (94-21)
    tmp = importdata("data\L2\dvchla_94-21_200.txt");
    ax = L2_helper(tmp,mld,dcm,50,testSel,"ad",[-60 60],[7 19]);
    sgtitle("HPLC Divinyl chlorophyll a (94-21): L2"+tmpT);
    exportgraphics(ax,"figures/L2/bottle/log/dvchla" + tmpT + ".png");
    % save("output\L2\dvchla.mat","p","ks","obs","sk","ku");
    clearvars -except mld dcm tmpT testSel;
    
    % HPLC Chlorophyll b (88-21)
    tmp = importdata("data\L2\chlb_88-21_200.txt");
    ax = L2_helper(tmp,mld,dcm,50,testSel,"ad",[-60 60],[7 19]);
    sgtitle("L2");
    exportgraphics(ax,"figures/L2/bottle/log/chlb" + tmpT + ".png");
    % save("output\L2\chlb.mat","p","ks","obs","sk","ku");
    clearvars -except mld dcm tmpT testSel;
    
    % HPLC Chlorophyll C1 + C2 + C3 (88-21)
    tmp = importdata("data\L2\chl123_88-21_200.txt");
    ax = L2_helper(tmp,mld,dcm,50,testSel,"ad",[-60 60],[7 19]);
    sgtitle("L2");
    exportgraphics(ax,"figures/L2/bottle/log/chl123" + tmpT + ".png");
    % save("output\L2\chl123.mat","p","ks","obs","sk","ku");
    clearvars -except mld dcm tmpT testSel;
    
    % HPLC alpha-Carotene (94-21)
    tmp = importdata("data\L2\acar_94-21_200.txt");
    ax = L2_helper(tmp,mld,dcm,50,testSel,"ad",[-60 60],[7 19]);
    sgtitle("HPLC alpha-Carotene (94-21): L2"+tmpT);
    exportgraphics(ax,"figures/L2/bottle/log/acar" + tmpT + ".png");
    % save("output\L2\acar.mat","p","ks","obs","sk","ku");
    clearvars -except mld dcm tmpT testSel;
    
    % But19 (88-21)
    tmp = importdata("data\L2\but19_88-21_200.txt");
    ax = L2_helper(tmp,mld,dcm,50,testSel,"ad",[-60 60],[7 19]);
    sgtitle("But19 (88-21): L2"+tmpT);
    exportgraphics(ax,"figures/L2/bottle/log/but19" + tmpT + ".png");
    % save("output\L2\but19.mat","p","ks","obs","sk","ku");
    clearvars -except mld dcm tmpT testSel;
    
    % HPLC 19 Hexanoyloxyfucoxanthin (88-21)
    tmp = importdata("data\L2\hex19_88-21_200.txt");
    ax = L2_helper(tmp,mld,dcm,50,testSel,"ad",[-60 60],[7 19]);
    sgtitle("HPLC 19' Hexanoyloxyfucoxanthin (88-21): L2"+tmpT);
    exportgraphics(ax,"figures/L2/bottle/log/hex19" + tmpT + ".png");
    % save("output\L2\hex19.mat","p","ks","obs","sk","ku");
    clearvars -except mld dcm tmpT testSel;
    
    % HPLC Zeaxanthin (88-21)
    tmp = importdata("data\L2\zeaxan_88-21_200.txt");
    ax = L2_helper(tmp,mld,dcm,50,testSel,"ad",[-60 60],[7 19]);
    sgtitle("HPLC Zeaxanthin (88-21): L2"+tmpT);
    exportgraphics(ax,"figures/L2/bottle/log/zeaxan" + tmpT + ".png");
    % save("output\L2\zeaxan.mat","p","ks","obs","sk","ku");
    clearvars -except mld dcm tmpT testSel;
    
    % Particulate Carbon (89-21)
    tmpx = "";
    tmp = importdata("data\L2\parc_89-21_200.txt");
    ax = L2_helper(tmp,mld,dcm,50,testSel,"ad",[-80 80],[6 22]);
    sgtitle("L2"+tmpx);
    exportgraphics(ax,"figures/L2/bottle/log/pc" + tmpT + ".png");
    % save("output\L2\pc.mat","p","ks","obs","sk","ku");
    clearvars -except mld dcm tmpT testSel;
    
    % Particulate Nitrogen (89-21)
    tmp = importdata("data\L2\parn_89-21_200.txt");
    ax = L2_helper(tmp,mld,dcm,50,testSel,"ad",[-60 60],[8 20]);
    sgtitle("Particulate Nitrogen 89-21: L2"+tmpT);
    exportgraphics(ax,"figures/L2/bottle/log/pn" + tmpT + ".png");
    % save("output\L2\pn.mat","p","ks","obs","sk","ku");
    clearvars -except mld dcm tmpT testSel;
    
    % Nitrate + Nitrite (88-21)
    tmp = importdata("data/L2/nit2_88-21_200.txt");
    ax = L2_helper(tmp,mld,dcm,50,testSel,"ad",[-60 60],[7 19]);
    sgtitle("Nitrate + Nitrite 88-21: L2"+tmpT);
    exportgraphics(ax,"figures/L2/bottle/log/nit2" + tmpT + ".png");
    % save("output\L2\nit2.mat","p","ks","obs","sk","ku");
    clearvars -except mld dcm tmpT testSel;
    
    % Phosphate (88-21)
    tmp = importdata("data/L2/pho_88-21_200.txt");
    ax = L2_helper(tmp,mld,dcm,50,testSel,"ad",[-60 60],[7 19]);
    sgtitle("Phosphate 88-21: L2"+tmpT);
    exportgraphics(ax,"figures/L2/bottle/log/phos" + tmpT + ".png");
    % save("output\L2\phos.mat","p","ks","obs","sk","ku");
    clearvars -except mld dcm tmpT testSel;
    
    % Oxygen (88-21)
    tmp = importdata("data/L2/oxy_88-21_200.txt");
    ax = L2_helper(tmp,mld,dcm,50,testSel,"ad",[-60 60],[6 18]);
    sgtitle("Dissolved Oxygen 88-21: L2"+tmpT);
    exportgraphics(ax,"figures/L2/bottle/log/boxy" + tmpT + ".png");
    % save("output\L2\boxy.mat","p","ks","obs","sk","ku");
    clearvars -except mld dcm tmpT testSel;
    
    % PProd Light-12 (89-22)
    tmp = importdata("data\L2\l12_89-22_200.txt");
    ax = L2_helper(tmp,mld,dcm,50,testSel,"ad",[-60 60],[6 18]);
    sgtitle("PProd Light-12 (89-22): L2"+tmpT);
    exportgraphics(ax,"figures/L2/bottle/log/l12" + tmpT + ".png");
    % save("output\L2\l12.mat","p","ks","obs","sk","ku");
    clearvars -except mld dcm tmpT testSel;
end

%% K-S Unused

% % Prochlorococcus (05-21)
% tmp = importdata("data/L2/pro_05-21_200.txt");
% [ax,p,ks,obs,sk,ku] = L2_helper(tmp,mld,dcm,[2 22]);
% sgtitle("Prochlorococcus 05-21: L2");
% exportgraphics(ax,"figures/L2/bottle/log/notUsed/pbact21" + tmpT + ".png"); clear ax;
% save("output\L2\pbact21.mat","p","ks","obs","sk","ku");
% clear ax p ks obs sk ku;

% % Prochlorococcus (90-05)
% tmp = importdata("data/L2/pro_90-05_200.txt");
% [ax,p,ks,obs,sk,ku] = L2_helper(tmp,mld,dcm,[3 23]);
% sgtitle("Prochlorococcus 90-05: L2");
% exportgraphics(ax,"figures/L2/pbact05.png");
% save("output\L2\pbact05.mat","p","ks","obs","sk","ku");
% clearvars -except mld dcm tmpT;

% % Dissolved Inorganic Carbon (88-21)
% tmp = importdata("data/L2/dic_88-21_200.txt");
% [ax,p,ks,obs,sk,ku] = L2_helper(tmp,mld,dcm,[2 22]);
% sgtitle("Dissolved Inorganic Carbon 88-21: L2");
% exportgraphics(ax,"figures/L2/bottle/log/notUsed/dic" + tmpT + ".png");
% save("output\L2\dic.mat","p","ks","obs","sk","ku");
% clearvars -except mld dcm tmpT;

% % pH (92-21)
% tmp = importdata("data/L2/pH_92-21_200.txt");
% [ax,p,ks,obs,sk,ku] = L2_helper(tmp,mld,dcm,[2 22]);
% sgtitle("pH 92-21: L2");
% exportgraphics(ax,"figures/L2/bottle/log/notUsed/pH" + tmpT + ".png");
% save("output\L2\pH.mat","p","ks","obs","sk","ku");
% clearvars -except mld dcm tmpT;

% % Alkalinity (89-21)
% tmp = importdata("data/L2/alk_89-21_200.txt");
% [ax,p,ks,obs,sk,ku] = L2_helper(tmp,mld,dcm,[2 22]);
% sgtitle("Alkalinity 89-21: L2");
% exportgraphics(ax,"figures/L2/bottle/log/notUsed/alk" + tmpT + ".png");
% save("output\L2\alk.mat","p","ks","obs","sk","ku");
% clearvars -except mld dcm tmpT;

% % Silicate (88-22)
% tmp = importdata("data/L2/sil_88-22_200.txt");
% [ax,p,ks,obs,sk,ku] = L2_helper(tmp,mld,dcm,[3 23]);
% sgtitle("Silicate 88-22: L2");
% exportgraphics(ax,"figures/L2/bottle/log/notUsed/sil" + tmpT + ".png");
% save("output\L2\sil.mat","p","ks","obs","sk","ku");
% clearvars -except mld dcm tmpT;

% % Dissolved Organic Phosphorus (DOP) (88-01)
% tmp = importdata("data\L2\dop_88-01_200.txt");
% [ax,p,ks,obs,sk,ku] = L2_helper(tmp,mld,dcm);
% sgtitle("DOP 88-01: L2");
% exportgraphics(ax,"figures/L2/bottle/log/notUsed/dop" + tmpT + ".png");
% save("output\L2\dop.mat","p","ks","obs","sk","ku");
% clearvars -except mld dcm tmpT;

% % DON (88-17)
% tmp = importdata("data\L2\don_88-17_200.txt");
% [ax,p,ks,obs,sk,ku] = L2_helper(tmp,mld,dcm);
% sgtitle("DON 88-17: L2");
% exportgraphics(ax,"figures/L2/bottle/log/notUsed/don" + tmpT + ".png"); clear ax;
% save("output\L2\don.mat","p","ks","obs","sk","ku");
% clear ax p ks obs sk ku tmpT;

% % DOC (93-17)
% tmp = importdata("data\L2\doc_93-17_200.txt");
% [ax,p,ks,obs,sk,ku] = L2_helper(tmp,mld,dcm,[3 23]);
% sgtitle("DOC 88-17: L2");
% exportgraphics(ax,"figures/L2/bottle/log/notUsed/doc" + tmpT + ".png");
% save("output\L2\doc.mat","p","ks","obs","sk","ku");
% clearvars -except mld dcm tmpT;

% % TDP (88-01)
% tmp = importdata("data\L2\tdp_88-01_200.txt");
% [ax,p,ks,obs,sk,ku] = L2_helper(tmp,mld,dcm);
% sgtitle("TDP 88-01: L2");
% exportgraphics(ax,"figures/L2/bottle/log/notUsed/tdp" + tmpT + ".png");
% save("output\L2\tdp.mat","p","ks","obs","sk","ku");
% clearvars -except mld dcm tmpT;

% % Total Dissolved Nitrogen (TDN) (88-17)
% tmp = importdata("data\L2\tdn_88-17_200.txt");
% [ax,p,ks,obs,sk,ku] = L2_helper(tmp,mld,dcm);
% sgtitle("TDN 88-17: L2");
% exportgraphics(ax,"figures/L2/bottle/log/notUsed/tdn" + tmpT + ".png");
% save("output\L2\tdn.mat","p","ks","obs","sk","ku");
% clearvars -except mld dcm tmpT;

% % Particulate Phosphorus (11-21)
% tmp = importdata("data\L2\parp_11-21_200.txt");
% [ax,p,ks,obs,sk,ku] = L2_helper(tmp,mld,dcm);
% sgtitle("Particulate Phosphorus 11-21: L2");
% exportgraphics(ax,"figures/L2/pp" + tmpT + ".png");
% save("output\L2\pp.mat","p","ks","obs","sk","ku");
% clearvars -except mld dcm tmpT;

% % Particulate Silica (96-21)
% tmp = importdata("data\L2\pars_96-21_200.txt");
% [ax,p,ks,obs,sk,ku] = L2_helper(tmp,mld,dcm,[3 23]);
% sgtitle("Particulate Silica 96-21: L2");
% exportgraphics(ax,"figures/L2/bottle/log/notUsed/ps" + tmpT + ".png");
% save("output\L2\ps.mat","p","ks","obs","sk","ku");
% clearvars -except mld dcm tmpT;

% % Low-level Phosphorus (88-22)
% tmp = importdata("data\L2\llp_88-22_200.txt");
% [ax,p,ks,obs,sk,ku] = L2_helper(tmp,mld,dcm,[11 31]);
% sgtitle("Low-level Phosphorus 88-22: L2");
% exportgraphics(ax,"figures/L2/bottle/log/notUsed/llp" + tmpT + ".png");
% save("output\L2\llp.mat","p","ks","obs","sk","ku");
% clearvars -except mld dcm tmpT;

% % Low-level Nitrogen (89-22)
% tmp = importdata("data\L2\lln_89-22_200.txt");
% [ax,p,ks,obs,sk,ku] = L2_helper(tmp,mld,dcm,[11 31]);
% sgtitle("Low-level Nitrogen 89-22: L2");
% exportgraphics(ax,"figures/L2/bottle/log/notUsed/lln" + tmpT + ".png");
% save("output\L2\lln.mat","p","ks","obs","sk","ku"); clearvars -except mld dcm tmpT;

% % Fluorometric Chlorophyll (88-22)
% tmp = importdata("data\L2\chlFlu_88-22_200.txt");
% [ax,p,ks,obs,sk,ku] = L2_helper(tmp,mld,dcm,[11 31]);
% sgtitle("Fluorometric Chlorophyll 88-22: L2");
% exportgraphics(ax,"figures/L2/bottle/log/notUsed/chlaFluo" + tmpT + ".png");
% save("output\L2\chlaFluo.mat","p","ks","obs","sk","ku");
% clearvars -except mld dcm tmpT;

% % Phaeopigments (88-22)
% tmp = importdata("data\L2\pheo_88-22_200.txt");
% [ax,p,ks,obs,sk,ku] = L2_helper(tmp,mld,dcm,[11 31]);
% sgtitle("Phaeopigments 88-22: L2");
% exportgraphics(ax,"figures/L2/bottle/log/notUsed/phaeo" + tmpT + ".png");
% save("output\L2\phaeo.mat","p","ks","obs","sk","ku");
% clearvars -except mld dcm tmpT;

% % HPLC Chlorophyll C3 (88-21)
% tmp = importdata("data\L2\chl3_88-21_200.txt");
% [ax,p,ks,obs,sk,ku] = L2_helper(tmp,mld,dcm,[1 21]);
% sgtitle("Chlorophyll C3 (88-21): L2");
% exportgraphics(ax,"figures/L2/bottle/log/notUsed/chl3" + tmpT + ".png");
% save("output\L2\chl3.mat","p","ks","obs","sk","ku");
% clearvars -except mld dcm tmpT;

% % HPLC Chlorophyll C1 + C2 (88-21)
% tmp = importdata("data\L2\chl12_88-21_200.txt");
% [ax,p,ks,obs,sk,ku] = L2_helper(tmp,mld,dcm,[1 21]);
% sgtitle("Chlorophyll C1 + C2 (88-21): L2");
% exportgraphics(ax,"figures/L2/bottle/log/notUsed/chl12" + tmpT + ".png");
% save("output\L2\chl12.mat","p","ks","obs","sk","ku");
% clearvars -except mld dcm tmpT;

% % HPLC Fucoxanthin (88-21)
% tmp = importdata("data\L2\fuco_88-21_200.txt");
% [ax,p,ks,obs,sk,ku] = L2_helper(tmp,mld,dcm,[3 23]);
% sgtitle("HPLC Fucoxanthin (88-21): L2");
% exportgraphics(ax,"figures/L2/bottle/log/notUsed/fuco" + tmpT + ".png");
% save("output\L2\fuco.mat","p","ks","obs","sk","ku");
% clearvars -except mld dcm tmpT;

% % HPLC Prasinoxanthin (88-21)
% tmp = importdata("data\L2\prasino_88-21_200.txt");
% [ax,p,ks,obs,sk,ku] = L2_helper(tmp,mld,dcm,[3 23]);
% sgtitle("HPLC Prasinoxanthin (88-21): L2");
% exportgraphics(ax,"figures/L2/bottle/log/notUsed/prasino" + tmpT + ".png");
% save("output\L2\prasino.mat","p","ks","obs","sk","ku");
% clearvars -except mld dcm tmpT;

% % HPLC Diadinoxanthin (88-21)
% tmp = importdata("data\L2\diadino_88-21_200.txt");
% [ax,p,ks,obs,sk,ku] = L2_helper(tmp,mld,dcm,[3 23]);
% sgtitle("HPLC Diadinoxanthin (88-21): L2");
% exportgraphics(ax,"figures/L2/bottle/log/notUsed/diadino" + tmpT + ".png");
% save("output\L2\diadino.mat","p","ks","obs","sk","ku");
% clearvars -except mld dcm tmpT;

% % HPLC beta-Carotene (94-21)
% tmp = importdata("data\L2\bcar_94-21_200.txt");
% [ax,p,ks,obs,sk,ku] = L2_helper(tmp,mld,dcm,[3 23]);
% sgtitle("HPLC beta-Carotene (94-21): L2");
% exportgraphics(ax,"figures/L2/bottle/log/notUsed/bcar" + tmpT + ".png");
% save("output\L2\bcar.mat","p","ks","obs","sk","ku");
% clearvars -except mld dcm tmpT;

% % HPLC Carotenes (88-21)
% tmp = importdata("data\L2\caroten_88-21_200.txt");
% [ax,p,ks,obs,sk,ku] = L2_helper(tmp,mld,dcm,[3 23]);
% sgtitle("HPLC Carotenes (88-21): L2");
% exportgraphics(ax,"figures/L2/bottle/log/notUsed/caroten" + tmpT + ".png");
% save("output\L2\caroten.mat","p","ks","obs","sk","ku");
% clearvars -except mld dcm tmpT;

% % HPLC chlorophyllide a (94-21)
% tmp = importdata("data\L2\chlda_94-21_200.txt");
% [ax,p,ks,obs,sk,ku] = L2_helper(tmp,mld,dcm,[3 23]);
% sgtitle("HPLC chlorophyllide a (94-21): L2");
% exportgraphics(ax,"figures/L2/bottle/log/notUsed/chlda" + tmpT + ".png");
% save("output\L2\chlda.mat","p","ks","obs","sk","ku");
% clearvars -except mld dcm tmpT;

% % HPLC Violaxanthin (94-21)
% tmp = importdata("data\L2\viol_94-21_200.txt");
% [ax,p,ks,obs,sk,ku] = L2_helper(tmp,mld,dcm,[3 23]);
% sgtitle("HPLC Violaxanthin (94-21): L2");
% exportgraphics(ax,"figures/L2/bottle/log/notUsed/viol" + tmpT + ".png");
% save("output\L2\viol.mat","p","ks","obs","sk","ku");
% clearvars -except mld dcm tmpT;

% HPLC Lutein (94-21)
% tmp = importdata("data\L2\lutein_94-21_200.txt");
% [ax,p,ks,obs,sk,ku] = L2_helper(tmp,mld,dcm);
% sgtitle("HPLC Lutein (94-21): L2");
% exportgraphics(ax,"figures/L2/lutein.png");
% save("output\L2\lutein.mat","p","ks","obs","sk","ku");
% clearvars -except mld dcm tmpT;

% % Phycoerythrin 0.4u fraction (00-08)
% tmp = importdata("data\L2\pe4_00-08_200.txt");
% [ax,p,ks,obs,sk,ku] = L2_helper(tmp,mld,dcm);
% sgtitle("Phycoerythrin 0.4u fraction (00-08): L2");
% exportgraphics(ax,"figures/L2/pe4.png");
% save("output\L2\pe4.mat","p","ks","obs","sk","ku");
% clearvars -except mld dcm tmpT;

% % Phycoerythrin 5u fraction (00-08)
% tmp = importdata("data\L2\pe5_00-08_200.txt");
% [ax,p,ks,obs,sk,ku] = L2_helper(tmp,mld,dcm);
% sgtitle("Phycoerythrin 5u fraction (00-08): L2");
% exportgraphics(ax,"figures/L2/pe5.png");
% save("output\L2\pe5.mat","p","ks","obs","sk","ku");
% clearvars -except mld dcm tmpT;

% % Phycoerythrin 10u fraction (00-08)
% tmp = importdata("data\L2\pe10_00-08_200.txt");
% [ax,p,ks,obs,sk,ku] = L2_helper(tmp,mld,dcm);
% sgtitle("Phycoerythrin 10u fraction (00-08): L2");
% exportgraphics(ax,"figures/L2/pe10.png");
% save("output\L2\pe10.mat","p","ks","obs","sk","ku");
% clearvars -except mld dcm tmpT;

% % Heterotrophic Bacteria (05-21)
% tmp = importdata("data\L2\hbact_05-21_200.txt");
% [ax,p,ks,obs,sk,ku] = L2_helper(tmp,mld,dcm);
% sgtitle("Heterotrophic Bacteria (05-21): L2");
% exportgraphics(ax,"figures/L2/bottle/log/notUsed/hbact" + tmpT + ".png");
% save("output\L2\hbact.mat","p","ks","obs","sk","ku");
% clearvars -except mld dcm tmpT;

% % Prochlorococcus (05-21)
% tmp = importdata("data\L2\pbact_05-21_200.txt");
% [ax,p,ks,obs,sk,ku] = L2_helper(tmp,mld,dcm);
% sgtitle("Prochlorococcus (05-21): L2");
% exportgraphics(ax,"figures/L2/pbact.png");
% save("output\L2\pbact.mat","p","ks","obs","sk","ku");
% clearvars -except mld dcm tmpT;

% % Synechococcus (05-21)
% tmp = importdata("data\L2\sbact_05-21_200.txt");
% [ax,p,ks,obs,sk,ku] = L2_helper(tmp,mld,dcm);
% sgtitle("Synechococcus (05-21): L2");
% exportgraphics(ax,"figures/L2/sbact.png");
% save("output\L2\sbact.mat","p","ks","obs","sk","ku");
% clearvars -except mld dcm tmpT;

% % Picoeukaryotes (05-21)
% tmp = importdata("data\L2\ebact_05-21_200.txt");
% [ax,p,ks,obs,sk,ku] = L2_helper(tmp,mld,dcm);
% sgtitle("Picoeukaryotes (05-21): L2");
% exportgraphics(ax,"figures/L2/bottle/log/notUsed/ebact" + tmpT + ".png");
% save("output\L2\ebact.mat","p","ks","obs","sk","ku");
% clearvars -except mld dcm tmpT;

% % ATP (88-22)
% tmp = importdata("data\L2\atp_88-22_200.txt");
% [ax,p,ks,obs,sk,ku] = L2_helper(tmp,mld,dcm,[3 23]);
% sgtitle("ATP (88-22): L2");
% exportgraphics(ax,"figures/L2/bottle/log/notUsed/atp" + tmpT + ".png");
% save("output\L2\atp.mat","p","ks","obs","sk","ku");
% clearvars -except mld dcm tmpT;

% % Nitrous Oxide (93-01)
% tmp = importdata("data\L2\n2o_93-01_200.txt");
% [ax,p,ks,obs,sk,ku] = L2_helper(tmp,mld,dcm,[4 14]);
% sgtitle("Nitrous Oxide (93-01): L2");
% exportgraphics(ax,"figures/L2/n2o.png");
% save("output\L2\n2o.mat","p","ks","obs","sk","ku");
% clearvars -except mld dcm tmpT;

% % PProd Dark-12 (89-00)
% tmp = importdata("data\L2\d12_89-00_200.txt");
% [ax,p,ks,obs,sk,ku] = L2_helper(tmp,mld,dcm,[2 12]);
% sgtitle("PProd Dark-12 (89-00): L2");
% exportgraphics(ax,"figures/L2/d12.png");
% save("output\L2\d12.mat","p","ks","obs","sk","ku");
% clearvars -except mld dcm tmpT;