function [mN,mL,pN,pL,mG,mW,pG,pW] = testSkewnessBiasHelper(s,runs,dispFig,Nps,Lps,Gps,Wps)
% INPUTS
% s = sample size
% OUTPUTS
% sn = skewness of normally-distributed random data
% sl = skewness of lognormally-distributed random data

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
sn = nan(runs,1);
for i = 1:runs
    %a = randn(s,1);
    a = normrnd(Nps(1),Nps(2),[s 1]);
    sn(i) = skewness(a);
end

% LOGNORMAL case
sl = nan(runs,1);
for i = 1:runs
    b = lognrnd(Lps(1),Lps(2),[s 1]);
    sl(i) = skewness(b);
end

if nargout > 5
    % GAMMA case
    sg = nan(runs,1);
    for i = 1:runs
        %[0.5 2 5 9],[0.5 1 1.5 2]
        c = gamrnd(Gps(1),Gps(2),[s 1]);
        sg(i) = skewness(c);
    end
    
    % WEIBULL case
    sw = nan(runs,1);
    for i = 1:runs
        %0.5 1 1.5 2],[0.5 1.75 3.5 5]
        d = wblrnd(Wps(1),Wps(2),[s 1]);
        sw(i) = skewness(d);
    end

    mG = mean(sg);
    mW = mean(sw);

    pG = [prctile(sg,16) prctile(sg,84)];
    pW = [prctile(sw,16) prctile(sw,84)];
end

if dispFig == true
    ax = figure;
    plot(sn,DisplayName='Random Normal');
    hold on
    plot(sl,DisplayName='Random Lognormal');
    hold off
    legend();
    title('Skewness of randomly-generated arrays',sprintf('Sample Size = %d',s));
    str = "t" + s + "_.png";
    exportgraphics(ax,'figures/skewBias/' + str); clear ax;
end

% Get MEAN and STD of runs of skewness
mN = mean(sn);
mL = mean(sl);
pN = [prctile(sn,16) prctile(sn,84)];
pL = [prctile(sl,16) prctile(sl,84)];

end