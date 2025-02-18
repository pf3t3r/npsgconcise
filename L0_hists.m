% Examine histograms of bottle data at a single, binned depth.

close all;clc;clear;

%% IMPORT data.
% ONLY need to change variables in this subsection.
tmp = importdata("data\L0\dvchla_94-21_200.txt").data;
NAME = "dvchla-55";
xlab = "divinyl chl $a$ [ng/l]";
binVal = 55;

%% SPECIFY plot dimensions.

set(groot, 'defaultFigureUnits', 'centimeters', 'defaultFigurePosition', [5 5 20 10]);
set(0,'defaultAxesFontSize',16);

%% EXTRACT data.

p = tmp(:,4); 
X = tmp(:,5);
n = length(p);

%% BIN data.

pB = discretize(p,0:10:200);
pX = nan(n,1);
for i = 1:n
    pX(i) = pB(i)*10 - 5;
end
clear i n p pB tmp;

%% PLOT.
tmpX = X(pX==binVal);
tmpX(tmpX<=0) = [];

figure;
histfit(tmpX,40,"lognormal");
title(NAME); xlabel(xlab,Interpreter="latex");
legend("data","lognormal");
exportgraphics(gca,"figures/L0/bot/hist/" + NAME + ".png")

%% HYPOTHESIS test: A-D, Lillie.
[h,p] = adtest(tmpX,"Distribution","logn");
[hl,pl] = lillietest(log(tmpX),"Distr","norm");
% [hc,pc] = jbtest(tmpX);