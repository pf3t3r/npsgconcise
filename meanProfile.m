% This script shows the mean profile of fluorescence.

clear;clc;close all;
set(groot, 'defaultFigureUnits', 'centimeters', 'defaultFigurePosition', [3 5 10 10]);

% Import Fluorescence data.
data = importdata('data/hots-chloropigment.txt').data;

p = data(:,3);  % pressure
f = data(:,4);  % fluorescence

f(f==-9) = nan; % nan flag

% reshape the data into grids
pp = reshape(p,[101 329]);
ff = reshape(f,[101 329]);

% Find the mean and 5th/95th percentile of the fluorescence time series
ffm = mean(ff,2,"omitnan");
f5 = prctile(ff',5);
f95 = prctile(ff',95);

% Toggle show y-label and title (for paper we don't need either)
displayYLabelAndTitle = false;

% Plot the mean profile of fluorescence with the mean and confidence
% interval.
ax = figure;
plot(f,p,'.',Color=[0.8 0.8 0.8],DisplayName="raw data");
hold on
plot(ffm,pp(:,1),'-',"Color",[0 0 0],DisplayName="mean");
plot(f5,pp(:,1),'-',"Color",[0.5 0.5 0.5],DisplayName="5%");
plot(f95,pp(:,1),'-',"Color",[0.5 0.5 0.5],DisplayName="95%");
hold off
set(gca,"YDir","reverse");
legend();
xlabel("chl-$a$ fluorescence [$\mu$g/L]",Interpreter="latex");
if displayYLabelAndTitle == true
    title("L0 Fluorescence 1988-2022",Interpreter="latex");
    ylabel("Pressure [dbar]",Interpreter="latex");
    yticklabels({});
end
yticks(0:20:200);
exportgraphics(ax,'figures/L0/meanProfile.png');
