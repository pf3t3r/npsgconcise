% Examine histograms of CTD data at a single, binned depth.

close all;clc;clear;

%% IMPORT data.
% Change ONLY this section.

tmp = importdata("data\L0\chlaCtd_200.txt").data;   % import data
NAME = "chla-20";                                   % set name for file + plot title
xlab = "chloropigment [$\mu$g/l]";                        % set xlabel
pressure = 20;                                        % pressure at which to evaluate histogram
histBins = 45;                                          % no. of bins in histogram
showTitle = false;                                  % true = show figure title
showDist = true;                                    % true = show fitted distribution
multipleDist = true;                                % false = lognormal only;
                                                    % true = show logn./norm./weib./gamma

%% SPECIFY plot dimensions.

set(groot, 'defaultFigureUnits', 'centimeters', 'defaultFigurePosition', [5 5 18 10]);
set(0,'defaultAxesFontSize',16);

%% EXTRACT data.

p = tmp(:,3); 
X = tmp(:,4);

%% PLOT.
tmpX = X(p==pressure);
tmpX(tmpX<=0) = [];

figure;
if showDist == true
    if multipleDist == false
        h = histfit(tmpX,histBins,"lognormal");
        h(1).FaceColor = [0.8 0.8 0.8];
        legend("","lognormal");
    else
        h1 = histfit(tmpX,histBins,"lognormal");
        h1(1).FaceColor = [0.8 0.8 0.8]; h1(2).Color = "#1f78b4";
        hold on
        h2 = histfit(tmpX,histBins,"normal"); h2(2).Color = "#a6cee3";
        h3 = histfit(tmpX,histBins,"weibull"); h3(2).Color = "#b2df8a";
        h4 = histfit(tmpX,histBins,"gamma"); h4(2).Color = "#33a02c";
        hold off
        delete(h2(1)); delete(h3(1)); delete(h4(1));
        legend("","lognormal","normal","weibull","gamma");
    end
else
    h = histogram(tmpX,histBins);
    h(1).FaceColor = [0.8 0.8 0.8];
end
if showTitle == true
    title(NAME);
end
xlabel(xlab,Interpreter="latex");
xlim([0 inf]); % Constrain X to positive values only
exportgraphics(gca,"figures/L0/ctd/hist/" + NAME + ".png")

%% HYPOTHESIS test: A-D, Lillie.
[tmp1,tmp2] = adtest(tmpX,"Distribution","logn");
[Lil(1),Lil(2)] = lillietest(log(tmpX),"Distr","norm");
AD = [tmp1 tmp2];

clear h1 h2 h3 h4 h multipleDist showDist showTitle xlab NAME tmp1 tmp2;