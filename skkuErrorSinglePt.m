% Load skewness and kurtosis from a single-parameter file (e.g. chl-a). 
% Visualise this skewness and kurtosis at one depth.
% Compare this observational data with estimates of skewness and kurtosis
% from random data based on distributions using parameters estimated on
% that same data. (cf. inspectMomentsBias.m)

close all;clc;clear;

data = load("output\L1\chla.mat");
sk = data.Sk;
ku = data.Ku;
p = data.p;
obs = data.obs;

figure
plot(sk(3),ku(3),'*',DisplayName="Data",Color='k');
hold on
errorbar(0,3,0.3,0.3,0.2,0.2,'o','Color',[0.6509803921568628 0.807843137254902 0.8901960784313725],LineWidth=1.9,DisplayName="Normal");
errorbar(1,4.8,1.1,1.3,0.25,0.25,'o','Color','#1f78b4',LineWidth=1.9,DisplayName="Lognormal");
errorbar(0.7,3.7,0.7,0.5,0.7,0.5,'o','Color','#33a02c',LineWidth=1.9,DisplayName="Gamma");
errorbar(0.2,2.75,0.25,0.25,0.15,0.15,'o','Color','#b2df8a',LineWidth=1.9,DisplayName="Weibull");
hold off
legend();
xlabel("Skewness"); ylabel("Kurtosis");