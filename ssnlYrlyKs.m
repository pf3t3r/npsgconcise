% Calculate seasonal and yearly Kolmogorov-Smirnov statistics for bottle
% chl-a (L0).

clear; clc; close all;
addpath("baroneRoutines\");
set(groot,'defaultAxesXGrid','on');
set(groot,'defaultAxesYGrid','on');
set(groot, 'defaultFigureUnits', 'centimeters', 'defaultFigurePosition', [3 5 15 15]);
set(0,'defaultAxesFontSize',10);

%% Load Data
chloro = load("datafiles\chloro.mat"','chloro256').chloro256;
chloroL = load("datafiles\chloro.mat"','chloro256l').chloro256l;
time = load("datafiles\chloro.mat","time256").time256;

depthMeasurements = 129;
eulerianDepth = linspace(0,2*depthMeasurements,depthMeasurements);
lagrangianDepth = linspace(-128,128,depthMeasurements);

%% Seasonal KS
% Need to apply this night time -> done in nightTimeKs.m

% Give numerical ID to each month
timeByMonth = discretize(month(time),12);

% Find ID for each season in the timeseries
winter = [];
spring = [];
summer = [];
autumn = [];

% 12:329 = 1989 Oct -> 2021 Dec (complete - first year)
% 12:196 = 1989 Oct -> 2007 Dec (BB's PhD)
for i=12:329
    if timeByMonth(i) <= 3
        winter = [winter i];
    elseif timeByMonth(i) >= 4 && timeByMonth(i) <= 6
        spring = [spring i];
    elseif timeByMonth(i) >= 7 && timeByMonth(i) <= 9
        summer = [summer i];
    else
        autumn = [autumn i];
    end
end

%% Apply KS to seasonal data

% Try Eulerian first.
ksE_winter = zeros(5,depthMeasurements);
ksE_spring = zeros(5,depthMeasurements);
ksE_summer = zeros(5,depthMeasurements);
ksE_autumn = zeros(5,depthMeasurements);

% Try Lagrangian later.
ksL_winter = zeros(5,depthMeasurements);
ksL_spring = zeros(5,depthMeasurements);
ksL_summer = zeros(5,depthMeasurements);
ksL_autumn = zeros(5,depthMeasurements);

% Eulerian
for i = 1:depthMeasurements
    disp(i);
    tmp = chloro(i,winter);
    tmp(isnan(tmp) | tmp<=0) = [];
    [~,ksE_winter(:,i),~] = statsplot2(tmp,'noplot');
    tmp2 = chloro(i,spring);
    tmp2(isnan(tmp2) | tmp2<=0) = [];
    [~,ksE_spring(:,i),~] = statsplot2(tmp2,'noplot');
    tmp3 = chloro(i,summer);
    tmp3(isnan(tmp3) | tmp3<=0) = [];
    [~,ksE_summer(:,i),~] = statsplot2(tmp3,'noplot');
    tmp4 = chloro(i,autumn);
    tmp4(isnan(tmp4) | tmp4<=0) = [];
    [~,ksE_autumn(:,i),~] = statsplot2(tmp4,'noplot');
end
clear tmp;

% Lagrangian
for i = 1:depthMeasurements
    disp(i);
    tmp = chloroL(i,winter);
    tmp(isnan(tmp) | tmp<=0) = [];
    [~,ksL_winter(:,i),~] = statsplot2(tmp,'noplot');
    tmp2 = chloroL(i,spring);
    tmp2(isnan(tmp2) | tmp2<=0) = [];
    [~,ksL_spring(:,i),~] = statsplot2(tmp2,'noplot');
    tmp3 = chloroL(i,summer);
    tmp3(isnan(tmp3) | tmp3<=0) = [];
    [~,ksL_summer(:,i),~] = statsplot2(tmp3,'noplot');
    tmp4 = chloroL(i,autumn);
    tmp4(isnan(tmp4) | tmp4<=0) = [];
    [~,ksL_autumn(:,i),~] = statsplot2(tmp4,'noplot');
end

%% Seasonal KS: Eulerian Figures

set(groot, 'defaultFigureUnits', 'centimeters', 'defaultFigurePosition', [5 2 25 18]);

ax8 = figure;
% WINTER
subplot(1,4,1)
plot(ksE_winter(1,:),eulerianDepth,'o-','Color',[0 0 0],'DisplayName','Normal','LineWidth',1.2,'MarkerSize',4);
hold on
plot(ksE_winter(2,:),eulerianDepth,'+--','Color',[0 0 0],'DisplayName','Lognormal','LineWidth',1.2,'MarkerSize',4);
plot(ksE_winter(3,:),eulerianDepth,'xr-','DisplayName','Weibull','MarkerSize',4);
plot(ksE_winter(4,:),eulerianDepth,'r.--','LineStyle','--','DisplayName','Gamma','MarkerSize',4);
hold off
set(gca, 'YDir','reverse');
legend('Location','best');
ylim([0 250]);
xlabel('p-value');
ylabel('Pressure [db]');
title('Winter (JFM)');

% SPRING
subplot(1,4,2)
plot(ksE_spring(1,:),eulerianDepth,'o-','Color',[0 0 0],'DisplayName','Normal','LineWidth',1.2,'MarkerSize',4);
hold on
plot(ksE_spring(2,:),eulerianDepth,'+--','Color',[0 0 0],'DisplayName','Lognormal','LineWidth',1.2,'MarkerSize',4);
plot(ksE_spring(3,:),eulerianDepth,'xr-','DisplayName','Weibull','MarkerSize',4);
plot(ksE_spring(4,:),eulerianDepth,'r.--','LineStyle','--','DisplayName','Gamma','MarkerSize',4);
hold off
set(gca, 'YDir','reverse');
ylim([0 250]);
xlabel('p-value');
ylabel('Pressure [db]');
title('Spring (AMJ)');

% SUMMER
subplot(1,4,3)
plot(ksE_summer(1,:),eulerianDepth,'o-','Color',[0 0 0],'DisplayName','Normal','LineWidth',1.2,'MarkerSize',4);
hold on
plot(ksE_summer(2,:),eulerianDepth,'+--','Color',[0 0 0],'DisplayName','Lognormal','LineWidth',1.2,'MarkerSize',4);
plot(ksE_summer(3,:),eulerianDepth,'xr-','DisplayName','Weibull','MarkerSize',4);
plot(ksE_summer(4,:),eulerianDepth,'r.--','LineStyle','--','DisplayName','Gamma','MarkerSize',4);
hold off
set(gca, 'YDir','reverse');
ylim([0 250]);
xlabel('p-value');
ylabel('Pressure [db]');
title('Summer (JAS)');

% AUTUMN
subplot(1,4,4)
plot(ksE_autumn(1,:),eulerianDepth,'o-','Color',[0 0 0],'DisplayName','Normal','LineWidth',1.4,'MarkerSize',4);
hold on
plot(ksE_autumn(2,:),eulerianDepth,'+--','Color',[0 0 0],'DisplayName','Lognormal','LineWidth',1.4,'MarkerSize',4);
plot(ksE_autumn(3,:),eulerianDepth,'xr-','DisplayName','Weibull','MarkerSize',4);
plot(ksE_autumn(4,:),eulerianDepth,'r.--','LineStyle','--','DisplayName','Gamma','MarkerSize',4);
hold off
set(gca, 'YDir','reverse');
ylim([0 250]);
xlabel('p-value');
ylabel('Pressure [db]');
title('Autumn (OND)');

sgtitle('Kolmogorov-Smirnov Test: Seasonal, Eulerian (1989-2021)');
exportgraphics(ax8,'figures/ks_seasonal_eulerian_89-21.png');

%% Seasonal KS: Lagrangian Figures

ax9 = figure;
% WINTER
subplot(1,4,1)
plot(ksL_winter(1,:),lagrangianDepth,'o-','Color',[0 0 0],'DisplayName','Normal','LineWidth',1.4,'MarkerSize',4);
hold on
plot(ksL_winter(2,:),lagrangianDepth,'+--','Color',[0 0 0],'DisplayName','Lognormal','LineWidth',1.4,'MarkerSize',4);
plot(ksL_winter(3,:),lagrangianDepth,'xr-','DisplayName','Weibull','MarkerSize',4);
plot(ksL_winter(4,:),lagrangianDepth,'r.--','LineStyle','--','DisplayName','Gamma','MarkerSize',4);
hold off
legend('Location','best');
set(gca, 'YDir','reverse');
ylim([-120 120]);
xlabel('p-value');
ylabel('Depth [m]');
title('Winter (JFM)');

% SPRING
subplot(1,4,2)
plot(ksL_spring(1,:),lagrangianDepth,'o-','Color',[0 0 0],'DisplayName','Normal','LineWidth',1.4,'MarkerSize',4);
hold on
plot(ksL_spring(2,:),lagrangianDepth,'+--','Color',[0 0 0],'DisplayName','Lognormal','LineWidth',1.4,'MarkerSize',4);
plot(ksL_spring(3,:),lagrangianDepth,'xr-','DisplayName','Weibull','MarkerSize',4);
plot(ksL_spring(4,:),lagrangianDepth,'r.--','LineStyle','--','DisplayName','Gamma','MarkerSize',4);
hold 
set(gca, 'YDir','reverse');
ylim([-120 120]);
xlabel('p-value');
ylabel('Depth [m]');
title('Spring (AMJ)');

% SUMMER
subplot(1,4,3)
plot(ksL_summer(1,:),lagrangianDepth,'o-','Color',[0 0 0],'DisplayName','Normal','LineWidth',1.4,'MarkerSize',4);
hold on
plot(ksL_summer(2,:),lagrangianDepth,'+--','Color',[0 0 0],'DisplayName','Lognormal','LineWidth',1.4,'MarkerSize',4);
plot(ksL_summer(3,:),lagrangianDepth,'xr-','DisplayName','Weibull','MarkerSize',4);
plot(ksL_summer(4,:),lagrangianDepth,'r.--','LineStyle','--','DisplayName','Gamma','MarkerSize',4);
hold off
set(gca, 'YDir','reverse');
ylim([-120 120]);
xlabel('p-value');
ylabel('Depth [m]');
title('Summer (JAS)');

% AUTUMN
subplot(1,4,4)
plot(ksL_autumn(1,:),lagrangianDepth,'o-','Color',[0 0 0],'DisplayName','Normal','LineWidth',1.4,'MarkerSize',4);
hold on
plot(ksL_autumn(2,:),lagrangianDepth,'+--','Color',[0 0 0],'DisplayName','Lognormal','LineWidth',1.4,'MarkerSize',4);
plot(ksL_autumn(3,:),lagrangianDepth,'xr-','DisplayName','Weibull','MarkerSize',4);
plot(ksL_autumn(4,:),lagrangianDepth,'r.--','LineStyle','--','DisplayName','Gamma','MarkerSize',4);
hold off
set(gca, 'YDir','reverse');
ylim([-120 120]);
xlabel('p-value');
ylabel('Depth [m]');
title('Autumn (OND)');

sgtitle('Kolmogorov-Smirnov Test: Seasonal, Lagrangian (1989-2021)');
exportgraphics(ax9,'figures/ks_seasonal_lagrangian_89-21.png');

%% KS Test: stuff Benedetto just found!

% Eulerian
eulerianData = load('datafiles\hot_1773.mat');
fluoE = eulerianData.FLS_hot_1773(2:end,2:end);
pE = eulerianData.FLS_hot_1773(2:end,1);
ks_allcastE = nan(5,130);
for i = 1:130
    disp(i);
    tmp = fluoE(i,:);
    tmp(isnan(tmp) | tmp<=0) = [];
    if length(tmp) > 2 
        [~,ks_allcastE(:,i),~] = statsplot2(tmp,'noplot');
    end
end

% Lagrangian
lagrangianData = load('datafiles\lagr_1758.mat');
fluoL = lagrangianData.FLS_lagr_1758(2:end,2:end);
pL = lagrangianData.FLS_lagr_1758(2:end,1);
ks_allcastL = nan(5,221);
for i = 1:221
    disp(i);
    tmp = fluoL(i,:);
    tmp(isnan(tmp) | tmp<=0) = [];
    if length(tmp) > 2 
        [~,ks_allcastL(:,i),~] = statsplot2(tmp,'noplot');
    end
end

%% Plot Lagrangian of 1758 profiles

ax10 = figure;
subplot(1,2,1)
plot(ks_allcastE(1,:),pE,'o-','Color',[0 0 0],'DisplayName','Normal','LineWidth',1.4,'MarkerSize',4)
hold on
plot(ks_allcastE(2,:),pE,'+--','Color',[0 0 0],'DisplayName','Lognormal','LineWidth',1.4,'MarkerSize',4);
plot(ks_allcastE(3,:),pE,'xr-','DisplayName','Weibull','MarkerSize',4);
plot(ks_allcastE(4,:),pE,'r.--','LineStyle','--','DisplayName','Gamma','MarkerSize',4);
hold off
ylim([0 250]);
set(gca,'YDir','reverse');

title('Eulerian')

subplot(1,2,2)
plot(ks_allcastL(1,:),pL,'o-','Color',[0 0 0],'DisplayName','Normal','LineWidth',1.4,'MarkerSize',4)
hold on
plot(ks_allcastL(2,:),pL,'+--','Color',[0 0 0],'DisplayName','Lognormal','LineWidth',1.4,'MarkerSize',4);
plot(ks_allcastL(3,:),pL,'xr-','DisplayName','Weibull','MarkerSize',4);
plot(ks_allcastL(4,:),pL,'r.--','LineStyle','--','DisplayName','Gamma','MarkerSize',4);
hold off
ylim([-125 125]);
set(gca,'YDir','reverse');
title('Lagrangian');

ax10.Position = [3 3 20 15];
sgtitle('Yearly Kolmogorov-Smirnov Test');
exportgraphics(ax10,'figures/ks_89-07_allCasts.png');

%%
clear i axa axb axc axd ax1 ax2 ax3 ax4 ax5 ax6 ax7 ax8 ax9 ax10;