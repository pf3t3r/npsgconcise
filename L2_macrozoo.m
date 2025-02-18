% Script to output L2 results for macrozooplankton tows for statistical
% analysis.

clear; clc; close all;
addpath("baroneRoutines\"); addpath("func\");
set(groot, "defaultFigureUnits", "centimeters", "defaultFigurePosition", [3 3 28 15]);

% Load maximum mixed-layer depth 'MLD' and cruise-averaged deep chlorophyll
% maximum depth 'DCM'.
mld = load("mldVals.mat").maxMld; % single maximum per cruise
dcm = load("output/dcm.mat").dcm; % pDcm + sigmaDcm (all casts, crn 1-329)

% macrozooplankton
tmp = importdata('data/L0/macrozoo_94-22_200.txt');

crn = tmp.data(:,1);
p = tmp.data(:,5);
carbon = tmp.data(:,6);

%%
crn2 = crn(crn<=329);
p2 = p(crn<=329);
carbon2 = carbon(crn<=329);

stuff = nan(length(1:length(crn2(1):crn2(end))),1);

for i = 1:length(stuff)
    %length(unique(crn2))
    if crn2(i)
        disp(crn2(i));
        stuff(i) = min(p2(crn==i+51));
    end
end