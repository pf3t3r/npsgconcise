% Statistical Analysis of the mixed layer (L1) at Station ALOHA for chl-a, other pigments, and BGC variables.
% We import the data and run a hypothesis test on it with the respective
% null hypotheses of normal and lognormal. We use the Anderson-Darling test
% since this is both more powerful than similar tests such as
% Kolmogorov-Smirnov and more flexible than tests such as Shapiro-Wilks
% which did not easily allow for testing of other distributions.

clear; clc; close all; addpath("func\");
set(groot, "defaultFigureUnits", "centimeters", "defaultFigurePosition", [3 3 15 15]);

% Possible test cases.
principle = true;               % Main analysis: A-D Test complete time series
seasonal = false;               % Seasonal analysis: A-D test applied to 
                                % four seasons
fourdist = false;
showHistograms = false;
testSel = 2;                    % 2 = norm + logn (default), 
                                % 4 = norm + logn + weib + gamm
thresh = 30;                    % Threshold for A-D Test to be accepted
logAxes = true;                 % Output p-values as log values on x-axis
if logAxes == true
    lp = "log/";
else
    lp = "";
end

%% Calculate Mixed Layer Depth
% Specifically, find the Maximum Mixed Layer Depth recorded for each 
% cruise "maxMld".

ctdData = importdata("datafiles\ctd_iso_ALL.mat").ctd;
maxMld = nan(329,1);
for i = 1:329
    if ~isnan([ctdData(i).mld003])
        maxMld(i) = max([ctdData(i).mld003]);
    end
end

clear ctdData i;
save mldVals.mat maxMld;


%% Check histograms

if showHistograms == true
    tmp = importdata("data/L1/hplcChla_88-21_150.txt");
    
    [~,~,~,~,~,~,~,~,~,chla_ML,p_ML] = L1_helper(tmp,maxMld,thresh,testSel,"ad",logAxes,0,true);
    
    pB = round(p_ML,-1);   % bin the pressure
    
    % Plot.
    % For chla measured at pressures in range 5-15 dbar, use pB = 10
    % For range 15-25 dbar, use pB = 20, etc.
    figure;
    histogram(chla_ML(pB==10));
    
    % histfit variation
    figure
    histfit(chla_ML(pB==40),10,"lognormal");
end

%% Principal Analysis: A-D

if principle == true
    % A-D
    tmpT = "-ad";
   
    % HPLC chl-a
    tmp = importdata("data/L1/hplcChla_88-21_150.txt");
    ax = L1_helper(tmp,maxMld,thresh,testSel,"ad",logAxes);
    sgtitle("L1 chl-$a$","Interpreter","latex");
    exportgraphics(ax,"figures/L1/bottle/" + lp + "chla" + tmpT + ".png");
    
    % chl-b
    tmp = importdata("data/L1/chlb_88-21_150.txt");
    ax = L1_helper(tmp,maxMld,thresh,testSel,"ad",logAxes);
    sgtitle("L1 chl-$b$","Interpreter","latex");
    exportgraphics(ax,"figures/L1/bottle/" + lp + "chlb" + tmpT + ".png");
    
    % Particulate Carbon
    tmp = importdata("data/L1/parc_89-21_150.txt");
    ax = L1_helper(tmp,maxMld,thresh,testSel,"ad",logAxes);
    sgtitle("L1 Particulate Carbon","Interpreter","latex");
    exportgraphics(ax,"figures/L1/bottle/" + lp + "pc" + tmpT + ".png");
    clearvars -except maxMld lp thresh testSel seasonal fourdist;
end

%% Seasonal Analysis: A-D

if seasonal == true

    thresh = 30;
    set(groot, "defaultFigureUnits", "centimeters", "defaultFigurePosition", [3 3 10 15]);

    %%% WINTER
    tmpx = " Winter"; tmpT = "-ad-01";

    % chla
    tmp = importdata("data/L1/hplcChla_88-21_150.txt");
    ax = L1_helper(tmp,maxMld,thresh,testSel,"ad",true,1);
    sgtitle("L1 chl-$a$" + tmpx,"Interpreter","latex");
    exportgraphics(ax,"figures/L1/bottle/" + lp + "chla" + tmpT + ".png");
    clear tmp ax;
    
    % chlb (88-21)
    tmp = importdata("data\L1\chlb_88-21_150.txt");
    ax = L1_helper(tmp,maxMld,thresh,testSel,"ad",true,1);
    sgtitle("L1 chl-$b$"+tmpx,"Interpreter","latex");
    exportgraphics(ax,"figures/L1/bottle/" + lp + "chlb" + tmpT + ".png");
    clear tmp ax;
    
    % Particulate Carbon (89-21)
    tmp = importdata("data\L1\parc_89-21_150.txt");
    ax = L1_helper(tmp,maxMld,thresh,testSel,"ad",true,1);
    sgtitle("L1 Particulate Carbon"+tmpx,"Interpreter","latex");
    exportgraphics(ax,"figures/L1/bottle/" + lp + "pc" + tmpT + ".png");
    clear tmp ax;
    
    %%% SPRING
    tmpx = " Spring"; tmpT = "-ad-02";

    % chla
    tmp = importdata("data/L1/hplcChla_88-21_150.txt");
    ax = L1_helper(tmp,maxMld,thresh,testSel,"ad",true,2);
    sgtitle("L1 chl-$a$" + tmpx,"Interpreter","latex");
    exportgraphics(ax,"figures/L1/bottle/" + lp + "chla" + tmpT + ".png");
    clear tmp ax;
    
    % chlb (88-21)
    tmp = importdata("data\L1\chlb_88-21_150.txt");
    ax = L1_helper(tmp,maxMld,thresh,testSel,"ad",true,2);
    sgtitle("L1 chl-$b$"+tmpx,"Interpreter","latex");
    exportgraphics(ax,"figures/L1/bottle/" + lp + "chlb" + tmpT + ".png");
    clear tmp ax;
    
    % Particulate Carbon (89-21)
    tmp = importdata("data\L1\parc_89-21_150.txt");
    ax = L1_helper(tmp,maxMld,thresh,testSel,"ad",true,2);
    sgtitle("L1 Particulate Carbon"+tmpx,"Interpreter","latex");
    exportgraphics(ax,"figures/L1/bottle/" + lp + "pc" + tmpT + ".png");
    clear tmp ax;

    %%% SUMMER
    tmpx = " Summer"; tmpT = "-ad-03";

    % chla
    tmp = importdata("data/L1/hplcChla_88-21_150.txt");
    ax = L1_helper(tmp,maxMld,thresh,testSel,"ad",true,3);
    sgtitle("L1 chl-$a$" + tmpx,"Interpreter","latex");
    exportgraphics(ax,"figures/L1/bottle/" + lp + "chla" + tmpT + ".png");
    clear tmp ax;

    % chlb (88-21)
    tmp = importdata("data\L1\chlb_88-21_150.txt");
    ax = L1_helper(tmp,maxMld,thresh,testSel,"ad",true,3);
    sgtitle("L1 chl-$b$"+tmpx,"Interpreter","latex");
    exportgraphics(ax,"figures/L1/bottle/" + lp + "chlb" + tmpT + ".png");
    clear tmp ax;
    
    % Particulate Carbon (89-21)
    tmp = importdata("data\L1\parc_89-21_150.txt");
    ax = L1_helper(tmp,maxMld,thresh,testSel,"ad",true,3);
    sgtitle("L1 Particulate Carbon"+tmpx,"Interpreter","latex");
    exportgraphics(ax,"figures/L1/bottle/" + lp + "pc" + tmpT + ".png");
    clear tmp ax;
    
    %%% AUTUMN
    tmpx = " Autumn"; tmpT = "-ad-04";

    % chla
    tmp = importdata("data/L1/hplcChla_88-21_150.txt");
    ax = L1_helper(tmp,maxMld,thresh,testSel,"ad",true,4);
    sgtitle("L1 chl-$a$" + tmpx,"Interpreter","latex");
    exportgraphics(ax,"figures/L1/bottle/" + lp + "chla" + tmpT + ".png");
    clear tmp ax;
    
    % chlb (88-21)
    tmp = importdata("data\L1\chlb_88-21_150.txt");
    ax = L1_helper(tmp,maxMld,thresh,testSel,"ad",true,4);
    sgtitle("L1 chl-$b$"+tmpx,"Interpreter","latex");
    exportgraphics(ax,"figures/L1/bottle/" + lp + "chlb" + tmpT + ".png");
    clear tmp ax;
    
    % Particulate Carbon (89-21)
    tmp = importdata("data\L1\parc_89-21_150.txt");
    ax = L1_helper(tmp,maxMld,thresh,testSel,"ad",true,4);
    sgtitle("L1 Particulate Carbon"+tmpx,"Interpreter","latex");
    exportgraphics(ax,"figures/L1/bottle/" + lp + "pc" + tmpT + ".png");
    clear tmp ax;
end

%% A-D but with Four Distributions
if fourdist == true
    % A-D
    tmpT = "-ad-4dist";
    thresh = 30;
    
    % HPLC chl-a
    tmp = importdata("data/L1/hplcChla_88-21_150.txt");
    ax = L1_helper(tmp,maxMld,thresh,4,"ad");
    sgtitle("L1 chl-$a$: four distributions","Interpreter","latex");
    exportgraphics(ax,"figures/L1/bottle/" + lp + "chla" + tmpT + ".png");
end