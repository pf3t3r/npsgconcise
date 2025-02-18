function [mN,mL,pN,pL,mG,mW,pG,pW] = testKurtosisBiasHelper(s,runs,dispFig,Nps,Lps,Gps,Wps)
% INPUTS
% s = sample size
% runs = no. of random arrays generated in for loop
% dispFig = set true if you want to see the figure, otherwise false
% OUTPUTS
% kn = kurtosis of normally-distributed random data
% kl = kurtosis of lognormally-distributed random data
% mN = mean kurtosis of normally-distributed random data
% mL = mean kurtosis of lognormally-distributed random data
% pN = lower and upper percentile of the normal kurtoses
% pL = lower and upper percentile of the lognormal kurtoses

if nargout < 5
    mG = []; mW = [];
    pG = []; pW = [];
end

if nargin < 4
    Nps = [1 0.5];
    Lps = [0 0.5];
    Gps = [2 0.5];
    Wps = [1 3];
end

if nargin < 3
    dispFig = false;
end

if nargin < 2
    runs = 100;
end

% NORMAL case
kn = nan(runs,1);
for i = 1:runs
    %a = randn(s,1);
    a = normrnd(Nps(1),Nps(2),[s 1]);
    kn(i) = kurtosis(a);
end

% LOGNORMAL case
kl = nan(runs,1);
for i = 1:runs
    b = lognrnd(Lps(1),Lps(2),[s 1]);
    kl(i) = kurtosis(b);
end

if nargout > 5
    % GAMMA case
    kg = nan(runs,1);
    for i = 1:runs
        %[0.5 2 5 9],[0.5 1 1.5 2]
        c = gamrnd(Gps(1),Gps(2),[s 1]);
        kg(i) = kurtosis(c);
    end
    
    % WEIBULL case
    kw = nan(runs,1);
    for i = 1:runs
        %0.5 1 1.5 2],[0.5 1.75 3.5 5]
        d = wblrnd(Wps(1),Wps(2),[s 1]);
        kw(i) = kurtosis(d);
    end

    mG = mean(kg);
    mW = mean(kw);

    pG = [prctile(kg,16) prctile(kg,84)];
    pW = [prctile(kw,16) prctile(kw,84)];
end

if dispFig == true
    ax = figure;
    plot(kn,DisplayName='Random Normal');
    hold on
    plot(kl,DisplayName='Random Lognormal');
    %plot(kg,DisplayName='Random Gamma');
    %plot(kw,DisplayName='Random Weibull');
    hold off
    legend();
    title('Kurtosis of randomly-generated arrays',sprintf('Sample Size = %d',s));
    str = "test" + s + "bias.png";
    exportgraphics(ax,'figures/kurtBias/' + str); clear ax;
end

% Get MEAN and PRCTILE of runs of kurtosis
mN = mean(kn);
mL = mean(kl);

pN = [prctile(kn,16) prctile(kn,84)];
pL = [prctile(kl,16) prctile(kl,84)];


end