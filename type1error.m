% This code evaluates how often a set of hypothesis tests reject the null
% hypothesis (that the data is lognormal) incorrectly.
% X0 - X5 represent the Type 1 Error Rate (% of incorrect H0 rejections).
% The higher this rate is the less reliable the test is.

close all;clc;clear;
addpath 'C:\Users\pfarrell\AppData\Roaming\MathWorks\MATLAB Add-Ons\Functions\Shapiro-Wilk and Shapiro-Francia normality tests'

%%

expRun = 2;
tmpMin = nan(expRun,1);
tmpMinId = nan(expRun,1);

for i = 1:expRun

    r = 5000;
    s = 250;
    
    for j = 1:r
        tmp = randn(s,1);
        [h(j),p(j)] = kstest(tmp);
        [h1(j),p1(j)] = lillietest(tmp);
        [h2(j),p2(j)] = adtest(tmp);
        [h3(j),p3(j)] = swtest(tmp);
        [h4(j),p4(j)] = chi2gof(tmp);
        [h5(j),p5(j)] = jbtest(tmp);
    end
    
    X0 = length(find(h(h==1)))/r;
    X1 = length(find(h1(h1==1)))/r;
    X2 = length(find(h2(h2==1)))/r;
    X3 = length(find(h3(h3==1)))/r;
    X4 = length(find(h4(h4==1)))/r;
    X5 = length(find(h5(h5==1)))/r;
    
    tmp1 = [X0 X1 X2 X3 X4 X5];
    [tmpMin(i),tmpMinId(i)] = min(tmp1);
    disp(i);
end
