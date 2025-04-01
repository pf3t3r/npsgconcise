% Script to output L2 bottle results for the statistical analysis.

clear; clc; close all;
addpath("baroneRoutines\"); addpath("func\");
set(groot, "defaultFigureUnits", "centimeters", "defaultFigurePosition", [3 3 15 15]);

% Load maximum mixed-layer depth 'MLD' and cruise-averaged deep chlorophyll
% maximum depth 'DCM'.
mld = load("mldVals.mat").maxMld; % single maximum per cruise
dcm = load("output/dcm.mat").dcm; % pDcm + sigmaDcm (all casts, crn 1-329)

% Possible test cases.
showHistograms = true;
principleAnalysis = false;  % main analysis
seasonalAnalysisKs = false;   % seasonality of statistics
seasonalAnalysisAd = false;   % seasonality of statistics
testSel = 2; % 2 = norm + logn; 4 = norm + logn + weib + gamm

% TEMPLATE (XX-YY)
% 1. Load data
% 2. Extract data beneath ML, centre around DCM
% 3. Calculate K-S (or A-D) p-value, Vuong LLR, skewness + kurtosis
% 4. Plot results

%% L2 Histograms
if showHistograms == true
    % Select dataset and bin for analysis
    D = "data/L2/hplcChla_88-21_200.txt";
    bin = 60;   % value here represents the midpoint of a bin,
                % e.g. 10 = bin from 5-15 dbar, 20 = bin from 15-25 dbar
    nameVar = "chl-$a$";

    % Import data; extract data within mixed layer.
    tmp = importdata(D);
    [~,X_L2,~,~,~,~,~,~,~,~,~,~,~,~,~,p_L2] = L2_helper(tmp,mld,dcm,30,2,"ad",[-60 60],[7 19],0,true);
    
    % Plot with histfit.
    figure
    histfit(X_L2(p_L2==bin),10,"lognormal"); title(""+nameVar+": "+bin+" dbar",Interpreter="latex");
    grid on
end
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

    
    % HPLC Chlorophyll b (88-21)
    tmp = importdata("data\L2\chlb_88-21_200.txt");
    ax = L2_helper(tmp,mld,dcm,thresh,testSel,"ks",[-60 60],[5 17],1);
    sgtitle("HPLC Chlorophyll b (88-21): L2"+tmpT);
    exportgraphics(ax,"figures/L2/bottle/log/chlb" + tmpT + ".png");
    clear tmp ax;
    
    % Particulate Carbon (89-21)
    tmp = importdata("data\L2\parc_89-21_200.txt");
    ax = L2_helper(tmp,mld,dcm,thresh,testSel,"ks",[-60 60],[2 14],1);
    sgtitle("Particulate Carbon 89-21: L2" + tmpT);
    exportgraphics(ax,"figures/L2/bottle/log/pc" + tmpT + ".png");
    clear tmp ax;
    

    %%% SPRING
    tmpT = "-02";

    % chla
    tmp = importdata("data/L2/hplcChla_88-21_200.txt");
    ax = L2_helper(tmp,mld,dcm,thresh,testSel,"ks",[-60 60],[7 19],2);
    sgtitle("[Chl a] 88-21: L2"+tmpT);
    exportgraphics(ax,"figures/L2/bottle/log/chla" + tmpT + ".png");
    clear tmp ax;

    % HPLC Chlorophyll b (88-21)
    tmp = importdata("data\L2\chlb_88-21_200.txt");
    ax = L2_helper(tmp,mld,dcm,thresh,testSel,"ks",[-60 60],[7 19],2);
    sgtitle("HPLC Chlorophyll b (88-21): L2"+tmpT);
    exportgraphics(ax,"figures/L2/bottle/log/chlb" + tmpT + ".png");
    clear tmp ax;
    
    
    % Particulate Carbon (89-21)
    tmp = importdata("data\L2\parc_89-21_200.txt");
    ax = L2_helper(tmp,mld,dcm,thresh,testSel,"ks",[-60 60],[8 20],2);
    sgtitle("Particulate Carbon 89-21: L2" + tmpT);
    exportgraphics(ax,"figures/L2/bottle/log/pc" + tmpT + ".png");
    clear tmp ax;
    


    %%% SUMMER
    tmpT = "-03";

    % chla
    tmp = importdata("data/L2/hplcChla_88-21_200.txt");
    ax = L2_helper(tmp,mld,dcm,15,testSel,"ks",[-80 80],[3 19],3);
    sgtitle("[Chl a] 88-21: L2"+tmpT);
    exportgraphics(ax,"figures/L2/bottle/log/chla" + tmpT + ".png");
    clear tmp ax;
 
    % HPLC Chlorophyll b (88-21)
    tmp = importdata("data\L2\chlb_88-21_200.txt");
    ax = L2_helper(tmp,mld,dcm,thresh,testSel,"ks",[-60 60],[5 17],3);
    sgtitle("HPLC Chlorophyll b (88-21): L2"+tmpT);
    exportgraphics(ax,"figures/L2/bottle/log/chlb" + tmpT + ".png");
    clear tmp ax;
    
    
    % Particulate Carbon (89-21)
    tmp = importdata("data\L2\parc_89-21_200.txt");
    ax = L2_helper(tmp,mld,dcm,thresh,testSel,"ks",[-60 60],[5 17],3);
    sgtitle("Particulate Carbon 89-21: L2" + tmpT);
    exportgraphics(ax,"figures/L2/bottle/log/pc" + tmpT + ".png");
    clear tmp ax;
    

    %%% AUTUMN
    tmpT = "-04";
    
    % chla
    tmp = importdata("data/L2/hplcChla_88-21_200.txt");
    ax = L2_helper(tmp,mld,dcm,15,testSel,"ks",[-50 50],[1 11],4);
    sgtitle("[Chl a] 88-21: L2"+tmpT);
    exportgraphics(ax,"figures/L2/bottle/log/chla" + tmpT + ".png");
    clear tmp ax;

    
    % HPLC Chlorophyll b (88-21)
    tmp = importdata("data\L2\chlb_88-21_200.txt");
    ax = L2_helper(tmp,mld,dcm,thresh,testSel,"ks",[-50 50],[1 11],4);
    sgtitle("HPLC Chlorophyll b (88-21): L2"+tmpT);
    exportgraphics(ax,"figures/L2/bottle/log/chlb" + tmpT + ".png");
    clear tmp ax;
    
    
    % Particulate Carbon (89-21)
    tmp = importdata("data\L2\parc_89-21_200.txt");
    ax = L2_helper(tmp,mld,dcm,thresh,testSel,"ks",[-50 50],[1 11],4);
    sgtitle("Particulate Carbon 89-21: L2" + tmpT);
    exportgraphics(ax,"figures/L2/bottle/log/pc" + tmpT + ".png");
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

    
    % Particulate Carbon (89-21)
    tmp = importdata("data\L2\parc_89-21_200.txt");
    ax = L2_helper(tmp,mld,dcm,thresh,testSel,"ad",[-60 60],[1 13],1);
    sgtitle("L2 Winter");
    exportgraphics(ax,"figures/L2/bottle/log/pc" + tmpT + ".png");
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

    
    
    % HPLC Chlorophyll b (88-21)
    tmp = importdata("data\L2\chlb_88-21_200.txt");
    ax = L2_helper(tmp,mld,dcm,thresh,testSel,"ad",[-60 60],[7 19],2);
    sgtitle("HPLC Chlorophyll b (88-21): L2"+tmpT);
    exportgraphics(ax,"figures/L2/bottle/log/chlb" + tmpT + ".png");
    clear tmp ax;
    
   
    % Particulate Carbon (89-21)
    tmp = importdata("data\L2\parc_89-21_200.txt");
    ax = L2_helper(tmp,mld,dcm,thresh,testSel,"ad",[-60 60],[8 20],2);
    sgtitle("L2 Spring");
    exportgraphics(ax,"figures/L2/bottle/log/pc" + tmpT + ".png");
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

   
    
    % HPLC Chlorophyll b (88-21)
    tmp = importdata("data\L2\chlb_88-21_200.txt");
    ax = L2_helper(tmp,mld,dcm,thresh,testSel,"ad",[-60 60],[5 17],3);
    sgtitle("HPLC Chlorophyll b (88-21): L2"+tmpT);
    exportgraphics(ax,"figures/L2/bottle/log/chlb" + tmpT + ".png");
    clear tmp ax;
    
    
    
    % Particulate Carbon (89-21)
    tmp = importdata("data\L2\parc_89-21_200.txt");
    ax = L2_helper(tmp,mld,dcm,thresh,testSel,"ad",[-60 60],[5 17],3);
    sgtitle("L2 Summer");
    exportgraphics(ax,"figures/L2/bottle/log/pc" + tmpT + ".png");
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

    
    
    % HPLC Chlorophyll b (88-21)
    tmp = importdata("data\L2\chlb_88-21_200.txt");
    ax = L2_helper(tmp,mld,dcm,thresh,testSel,"ad",[-50 50],[1 11],4);
    sgtitle("HPLC Chlorophyll b (88-21): L2"+tmpT);
    exportgraphics(ax,"figures/L2/bottle/log/chlb" + tmpT + ".png");
    clear tmp ax;
    
    % Particulate Carbon (89-21)
    tmp = importdata("data\L2\parc_89-21_200.txt");
    ax = L2_helper(tmp,mld,dcm,thresh,testSel,"ad",[-50 50],[1 11],4);
    sgtitle("L2 Autumn");
    exportgraphics(ax,"figures/L2/bottle/log/pc" + tmpT + ".png");
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
       
    % HPLC Chlorophyll b (88-21)
    tmp = importdata("data\L2\chlb_88-21_200.txt");
    [ax,p,ks,obs,sk,ku] = L2_helper(tmp,mld,dcm,50,testSel,"ks",[-60 60],[7 19]);
    sgtitle("HPLC Chlorophyll b (88-21): L2");
    exportgraphics(ax,"figures/L2/bottle/log/chlb" + tmpT + ".png");
    save("output\L2\chlb.mat","p","ks","obs","sk","ku");
    clearvars -except mld dcm tmpT testSel;
    
    
    % Particulate Carbon (89-21)
    tmp = importdata("data\L2\parc_89-21_200.txt");
    [ax,p,ks,obs,sk,ku] = L2_helper(tmp,mld,dcm,50,testSel,"ks",[-80 80],[6 22]);
    sgtitle("Particulate Carbon 89-21: L2");
    exportgraphics(ax,"figures/L2/bottle/log/pc" + tmpT + ".png");
    save("output\L2\pc.mat","p","ks","obs","sk","ku");
    clearvars -except mld dcm tmpT testSel;
    
    
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
       
    % HPLC Chlorophyll b (88-21)
    tmp = importdata("data\L2\chlb_88-21_200.txt");
    ax = L2_helper(tmp,mld,dcm,50,testSel,"ad",[-60 60],[7 19]);
    sgtitle("L2");
    exportgraphics(ax,"figures/L2/bottle/log/chlb" + tmpT + ".png");
    % save("output\L2\chlb.mat","p","ks","obs","sk","ku");
    clearvars -except mld dcm tmpT testSel;
    
    % Particulate Carbon (89-21)
    tmpx = "";
    tmp = importdata("data\L2\parc_89-21_200.txt");
    ax = L2_helper(tmp,mld,dcm,50,testSel,"ad",[-80 80],[6 22]);
    sgtitle("L2"+tmpx);
    exportgraphics(ax,"figures/L2/bottle/log/pc" + tmpT + ".png");
    % save("output\L2\pc.mat","p","ks","obs","sk","ku");
    clearvars -except mld dcm tmpT testSel;
   
end