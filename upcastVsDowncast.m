% Compare upcast and downcast from individual casts. Upcasts are generally
% more reliable; it appears that downcasts occasionally record excess
% fluorescence at surface.

close all; clc; clear;

% Set Figure Parameters
set(groot,'defaultAxesXGrid','on');
set(groot,'defaultAxesYGrid','on');
set(groot, 'defaultFigureUnits', 'centimeters', 'defaultFigurePosition', [1 3 31 15]);
set(0,'defaultAxesFontSize',12);

addpath("baroneRoutines\");

%%
d = nan(6,515,10);

dplace = load("data/upcasts/h161a0204.2up");
d(1,1:length(dplace(:,1)),1:length(dplace(1,:))) = dplace;
clear dplace;

dplace = load("data/upcasts/h172a0209.2up");
d(2,1:length(dplace(:,1)),1:length(dplace(1,:))) = dplace;
clear dplace;

dplace = load("data/upcasts/h174a0204.2up");
d(3,1:length(dplace(:,1)),1:length(dplace(1,:))) = dplace;
clear dplace;

dplace = load("data/upcasts/h189a0204.2up");
d(4,1:length(dplace(:,1)),1:length(dplace(1,:))) = dplace;
clear dplace;

dplace = load("data/upcasts/h190a0207.2up");
d(5,1:length(dplace(:,1)),1:length(dplace(1,:))) = dplace;
clear dplace;

dplace = load("data/upcasts/h190a0215.2up");
d(6,1:length(dplace(:,1)),1:length(dplace(1,:))) = dplace;
clear dplace;

d = d(:,4:end,:);
d(d==-9) = NaN;

%% 

for i = 1:length(d(:,1,1))
    disp(i);
    figure;

    subplot(1,5,1)
    plot(d(i,:,2),d(i,:,1));
    hold on
    plot(d(i,:,8),d(i,:,1));
    hold off
    xlim([21 27]);
    ylabel('pressure [db]');
    xlabel('temperature [C]');
    set(gca,'YDir','reverse');

    subplot(1,5,2)
    plot(d(i,:,3),d(i,:,1));
    hold on
    plot(d(i,:,9),d(i,:,1));
    hold off
    xlim([34.5 35.5]);
    xlabel('salinity [g/kg]');
    set(gca,'YDir','reverse');

    subplot(1,5,3)
    plot(d(i,:,4),d(i,:,1));
    hold on
    plot(d(i,:,5),d(i,:,1));
    hold off
    xlim([180 230]);
    xlabel('nitrate?');
    set(gca,'YDir','reverse');

    subplot(1,5,4)
    plot(d(i,1:500,7),d(i,1:500,1));
    ylim([0 180]);
    xlabel('fluorescence?');
    set(gca,'YDir','reverse');
    
    subplot(1,5,5)
    plot(d(i,1:500,10),d(i,1:500,1));
    set(gca,'YDir','reverse');
    xlabel('microturbulent fourier phosphorus');
end
clear i;

% column 1 = pressure
% column 2 = temperature
% column 3 = salinity
% column 4 = nitrate?
% column 5 = nitrate?
% column 6 = ? (-99)
% column 7 = possibly fluorescence?
% column 8 = temperature (conservative?)
% column 9 = salinity (absolute?)
% column 10 = ?


%%
fl = load("datafiles\flagged.mat").flagged;

ax = figure;
for i = 1:6
    subplot(1,6,i)
    plot(d(i,1:90,7),d(i,1:90,1),'DisplayName','upcast');
    hold on
    plot(fl(1:90,i),d(i,1:90,1),'DisplayName','downcast');
    hold off
    ylim([0 180]);
    if i==1
        legend();
        ylabel('Pressure [db]');
    end
    xlabel('Fluorescence [\mugL^{-1}]');
    set(gca,'YDir','reverse');
end
sgtitle('Elevated surface fluorescence: comparison of six up/down casts');
exportgraphics(ax,'figures/flagme.png');