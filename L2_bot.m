% Statistical Analysis of the DCM layer (L2)* at Station ALOHA for chl-a, other pigments, and BGC variables.
% We import the data and run a hypothesis test on it with the respective
% null hypotheses of normal and lognormal. We use the Anderson-Darling (A-D)
% test since this is both more powerful than similar tests such as
% Kolmogorov-Smirnov (K-S) and more flexible than tests such as Shapiro-Wilks
% which did not easily allow for testing of other distributions. However, I
% have included K-S as an alternative here for comparison.

% *The DCM layer is defined as the region beneath the mixed layer that is
% centred on the Deep Chlorophyll Maximum (DCM).

clear; clc; close all;
addpath("func\");
set(groot, "defaultFigureUnits", "centimeters", "defaultFigurePosition", [3 3 15 15]);

% Options and Test Cases.
thresh = 30;                    % Set threshold: run A-D (or K-S) test only
                                % when no. of measurements at that depth 
                                % exceed this value.
testSel = 4;                    % 2 = norm + logn;
                                % 4 = norm + logn + weib + gamm
showHistograms = false;         % Show histogram analysis (T/F)
principleAnalysisAd = true;     % Show main results using A-D.
seasonalAnalysisAd = false;     % Show seasonal results using A-D.
principleAnalysisKs = false;    % Show main results using K-S.
seasonalAnalysisKs = false;     % Show seasonal results using K-S.

% Load maximum mixed-layer depth 'MLD' and cruise-averaged deep chlorophyll
% maximum depth 'DCM'.
mld = load("output/mldVals.mat").maxMld; % single maximum per cruise
dcm = load("output/dcm.mat").dcm; % pDcm + sigmaDcm (all casts, crn 1-329)

%% Analyse histogams.
if showHistograms == true
    % Select dataset and bin for analysis
    D = "input/L2/hplcChla_88-21_200.txt";
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


%% Principal Analysis.
if principleAnalysisAd == true
    % A-D
    tmpX = ": L2";
    
    % Chlorophyll a (88-21)
    tmp = importdata("input/L2/hplcChla_88-21_200.txt");
    L2_helper(tmp,mld,dcm,30,testSel,"ad",[-60 60],[7 19]);
    sgtitle("chl-$a$"+tmpX,"Interpreter","latex");
       
    % HPLC Chlorophyll b (88-21)
    tmp = importdata("input\L2\chlb_88-21_200.txt");
    L2_helper(tmp,mld,dcm,50,testSel,"ad",[-60 60],[7 19]);
    sgtitle("chl-$b$"+tmpX,"Interpreter","latex");
        
    % Particulate Carbon (89-21)
    tmp = importdata("input\L2\parc_89-21_200.txt");
    L2_helper(tmp,mld,dcm,50,testSel,"ad",[-80 80],[6 22]);
    sgtitle("Particulate Carbon"+tmpX,"Interpreter","latex");
   
end


%% Seasonal Analysis
if seasonalAnalysisAd == true

    set(groot, "defaultFigureUnits", "centimeters", "defaultFigurePosition", [3 3 10 15]);

    %%% WINTER
    tmpX = ": L2 Winter";

    % Chlorophyll a (88-21)
    tmp = importdata("input/L2/hplcChla_88-21_200.txt");
    L2_helper(tmp,mld,dcm,thresh,testSel,"ad",[-50 50],[4 14],1);
    sgtitle("chl-$a$"+tmpX,"Interpreter","latex");
    clear tmp;

    % HPLC Chlorophyll b (88-21)
    tmp = importdata("input\L2\chlb_88-21_200.txt");
    L2_helper(tmp,mld,dcm,thresh,testSel,"ad",[-60 60],[3 15],1);
    sgtitle("chl-$b$"+tmpX,"Interpreter","latex");
    clear tmp;

    % Particulate Carbon (89-21)
    tmp = importdata("input\L2\parc_89-21_200.txt");
    L2_helper(tmp,mld,dcm,thresh,testSel,"ad",[-60 60],[1 13],1);
    sgtitle("Particulate Carbon"+tmpX,"Interpreter","latex");
    clear tmp;
    

    %%% SPRING
    tmpX = ": L2 Spring";

    % chla
    tmp = importdata("input/L2/hplcChla_88-21_200.txt");
    L2_helper(tmp,mld,dcm,thresh,testSel,"ad",[-50 50],[7 17],2);
    sgtitle("chl-$a$"+tmpX,"Interpreter","latex");
    clear tmp;

    % HPLC Chlorophyll b (88-21)
    tmp = importdata("input\L2\chlb_88-21_200.txt");
    L2_helper(tmp,mld,dcm,thresh,testSel,"ad",[-60 60],[6 18],2);
    sgtitle("chl-$b$"+tmpX,"Interpreter","latex");
    clear tmp;
   
    % Particulate Carbon (89-21)
    tmp = importdata("input\L2\parc_89-21_200.txt");
    L2_helper(tmp,mld,dcm,thresh,testSel,"ad",[-60 60],[8 20],2);
    sgtitle("Particulate Carbon"+tmpX,"Interpreter","latex");
    clear tmp;
    

    %%% SUMMER
    tmpX = ": L2 Summer";

    % chla
    tmp = importdata("input/L2/hplcChla_88-21_200.txt");
    L2_helper(tmp,mld,dcm,thresh,testSel,"ad",[-50 50],[6 16],3);
    sgtitle("chl-$a$"+tmpX,"Interpreter","latex");
    clear tmp;
    
    % HPLC Chlorophyll b (88-21)
    tmp = importdata("input\L2\chlb_88-21_200.txt");
    L2_helper(tmp,mld,dcm,thresh,testSel,"ad",[-60 60],[5 17],3);
    sgtitle("chl-$b$"+tmpX,"Interpreter","latex");
    clear tmp;
    
    % Particulate Carbon (89-21)
    tmp = importdata("input\L2\parc_89-21_200.txt");
    L2_helper(tmp,mld,dcm,thresh,testSel,"ad",[-60 60],[5 17],3);
    sgtitle("Particulate Carbon"+tmpX,"Interpreter","latex");
    clear tmp;
    

    %%% AUTUMN
    tmpX = ": L2 Autumn";
    
    % chla
    tmp = importdata("input/L2/hplcChla_88-21_200.txt");
    L2_helper(tmp,mld,dcm,thresh,testSel,"ad",[-50 50],[3 13],4);
    sgtitle("chl-$a$"+tmpX,"Interpreter","latex");
    clear tmp;
    
    % HPLC Chlorophyll b (88-21)
    tmp = importdata("input\L2\chlb_88-21_200.txt");
    L2_helper(tmp,mld,dcm,thresh,testSel,"ad",[-50 50],[3 13],4);
    sgtitle("chl-$b$"+tmpX,"Interpreter","latex");
    clear tmp;
    
    % Particulate Carbon (89-21)
    tmp = importdata("input\L2\parc_89-21_200.txt");
    L2_helper(tmp,mld,dcm,thresh,testSel,"ad",[-50 50],[5 15],4);
    sgtitle("Particulate Carbon"+tmpX,"Interpreter","latex");
    clear tmp;
    
end

%% Additional analysis using K-S.
if principleAnalysisKs == true
    tmpX = ": L2";
    
    % Chlorophyll a (88-21)
    tmp = importdata("input/L2/hplcChla_88-21_200.txt");
    L2_helper(tmp,mld,dcm,50,testSel,"ks",[-60 60],[7 19]);
    sgtitle("chl-$a$"+tmpX,"Interpreter","latex");

    % HPLC Chlorophyll b (88-21)
    tmp = importdata("input\L2\chlb_88-21_200.txt");
    L2_helper(tmp,mld,dcm,50,testSel,"ks",[-60 60],[7 19]);
    sgtitle("chl-$b$"+tmpX,"Interpreter","latex");
    
    % Particulate Carbon (89-21)
    tmp = importdata("input\L2\parc_89-21_200.txt");
    L2_helper(tmp,mld,dcm,50,testSel,"ks",[-80 80],[6 22]);
    sgtitle("Particulate Carbon"+tmpX,"Interpreter","latex");
end

if seasonalAnalysisKs == true
    
    %%% WINTER
    tmpX = ": L2 Winter";

    % Chlorophyll a (88-21)
    tmp = importdata("input/L2/hplcChla_88-21_200.txt");
    L2_helper(tmp,mld,dcm,thresh,testSel,"ks",[-60 60],[3 15],1);
    sgtitle("chl-$a$"+tmpX,"Interpreter","latex");
    clear tmp;
    
    % HPLC Chlorophyll b (88-21)
    tmp = importdata("input\L2\chlb_88-21_200.txt");
    L2_helper(tmp,mld,dcm,thresh,testSel,"ks",[-60 60],[3 15],1);
    sgtitle("chl-$b$"+tmpX,"Interpreter","latex");
    clear tmp;
    
    % Particulate Carbon (89-21)
    tmp = importdata("input\L2\parc_89-21_200.txt");
    ax = L2_helper(tmp,mld,dcm,thresh,testSel,"ks",[-60 60],[1 13],1);
    sgtitle("Particulate Carbon" + tmpX,"Interpreter","latex");
    clear tmp;
    

    %%% SPRING
    tmpX = ": L2 Spring";

    % chla
    tmp = importdata("input/L2/hplcChla_88-21_200.txt");
    L2_helper(tmp,mld,dcm,thresh,testSel,"ks",[-60 60],[6 18],2);
    sgtitle("chl-$a$"+tmpX,"Interpreter","latex");
    clear tmp;

    % HPLC Chlorophyll b (88-21)
    tmp = importdata("input\L2\chlb_88-21_200.txt");
    L2_helper(tmp,mld,dcm,thresh,testSel,"ks",[-60 60],[6 18],2);
    sgtitle("chl-$b$"+tmpX,"Interpreter","latex");
    clear tmp;
    
    % Particulate Carbon (89-21)
    tmp = importdata("input\L2\parc_89-21_200.txt");
    L2_helper(tmp,mld,dcm,thresh,testSel,"ks",[-60 60],[8 20],2);
    sgtitle("Particulate Carbon" + tmpX,"Interpreter","latex");
    clear tmp;
    

    %%% SUMMER
    tmpX = ": L2 Summer";

    % chla
    tmp = importdata("input/L2/hplcChla_88-21_200.txt");
    L2_helper(tmp,mld,dcm,15,testSel,"ks",[-80 80],[3 19],3);
    sgtitle("chl-$a$"+tmpX,"Interpreter","latex");
    clear tmp;
 
    % HPLC Chlorophyll b (88-21)
    tmp = importdata("input\L2\chlb_88-21_200.txt");
    L2_helper(tmp,mld,dcm,thresh,testSel,"ks",[-60 60],[5 17],3);
    sgtitle("chl-$b$"+tmpX,"Interpreter","latex");
    clear tmp;
    
    % Particulate Carbon (89-21)
    tmp = importdata("input\L2\parc_89-21_200.txt");
    L2_helper(tmp,mld,dcm,thresh,testSel,"ks",[-60 60],[5 17],3);
    sgtitle("Particulate Carbon" + tmpX,"Interpreter","latex");
    clear tmp;
    

    %%% AUTUMN
    tmpX = ": L2 Autumn";
    
    % chla
    tmp = importdata("input/L2/hplcChla_88-21_200.txt");
    L2_helper(tmp,mld,dcm,15,testSel,"ks",[-50 50],[3 13],4);
    sgtitle("chl-$a$"+tmpX,"Interpreter","latex");
    clear tmp;
    
    % HPLC Chlorophyll b (88-21)
    tmp = importdata("input\L2\chlb_88-21_200.txt");
    L2_helper(tmp,mld,dcm,thresh,testSel,"ks",[-50 50],[3 13],4);
    sgtitle("chl-$b$"+tmpX,"Interpreter","latex");
    clear tmp;
    
    % Particulate Carbon (89-21)
    tmp = importdata("input\L2\parc_89-21_200.txt");
    L2_helper(tmp,mld,dcm,thresh,testSel,"ks",[-50 50],[5 15],4);
    sgtitle("Particulate Carbon" + tmpX,"Interpreter","latex");
    clear tmp;
    
end
