function ax = L0_ctdHelper(X,hypTest,logAxis)
%L0_CTDHELPER 

if nargin <3
    logAxis = false;
end
if nargin <2
    hypTest = "ks";
end

n = 101;
ks = nan(5,n); ad = nan(2,n);
sk = nan(1,n); ku = nan(1,n);
% rVC = nan(10,n); pVC = nan(10,n);
obs = nan(1,n);
pCtd = 0:2:200;
alphaKs = 0.05;

% 4. Calculate KS p-value, skewness, kurtosis
for i = 1:n
    % find concentration X_i at binned pressure i
    X_i = X(i,:);
    % apply KS test to X_i
    % change limit below to >3 to fix error with picoeu -> may change other
    % results
    X_i(isnan(X_i)) = [];
    if length(X_i) > 3
        if strcmp(hypTest,"ks")
            [~,ks(:,i),~] = statsplot2(X_i,'noplot');
        else
            [~,ad(1,i)] = adtest(X_i,"Distribution","norm");
            [~,ad(2,i)] = adtest(X_i,"Distribution","logn");
        end
        [rVC(:,i),pVC(:,i)] = bbvuong(X_i);
        sk(i) = skewness(X_i);
        ku(i) = kurtosis(X_i);
    end
    obs(i) = length(X_i);
    %clear X_i;
end

%%
% Lognormal family: generate theoretical skewness and kurtosis
sigTh = linspace(0,1,1000);
for i = 1:length(sigTh)
    skLogn(i) = (exp(sigTh(i)^2) + 2)*(sqrt(exp(sigTh(i)^2) - 1));
    kuLogn(i) = exp(4*sigTh(i)^2) + 2*exp(3*sigTh(i)^2) + 3*exp(2*sigTh(i)^2) - 3;
end

% Negative Distributions
skLognN = -skLogn;
kuLognN = kuLogn;

kurtLimB = 10; skewLimA = 0; skewLimB = 2.5;
if max(ku) > 10 & min(sk) < 0
    kurtLimB = max(ku) + 1;
    skewLimA = min(sk) - 0.1;
    skewLimB = max(sk) + 0.1;
elseif max(ku) > 10
    kurtLimB = max(ku) + 1;
    skewLimB = max(sk) + 0.1;
elseif min(sk) < 0 
    skewLimA = min(sk) - 0.1;
elseif max(sk) > 2.5
    kurtLimB = max(ku) + 1;
    skewLimB = max(sk) + 0.1;
else 
    kurtLimB = 10;
    skewLimA = 0;
    skewLimB = 2.5;
end

%% K-S SK-KU Figure
set(groot, 'defaultFigureUnits', 'centimeters', 'defaultFigurePosition', [3 3 18 15]);

% pXX = 5:10:195;
ax = figure;
subplot(1,2,1)
xline(0.05,DisplayName="\alpha");
hold on
if strcmp(hypTest,"ks")
    plot(ks(2,:),pCtd,'+--','Color','#1f78b4',LineWidth=1.5,MarkerSize=5,HandleVisibility='off');
    xlabel('K-S $p$-value','Interpreter','latex',FontSize=15);
else
    plot(ad(1,:),pCtd,'o-','Color','#c51b7d',LineWidth=1.5,MarkerSize=5,DisplayName="normal");
    plot(ad(2,:),pCtd,'o-','Color','#4d9221',LineWidth=1.5,MarkerSize=5,DisplayName="lognormal");
    xlabel('A-D $p$-value','Interpreter','latex',FontSize=15);
end
if logAxis == true
    set(gca, 'XScale', 'log');
    xline(0.005,'--',HandleVisibility='off');
    xline(0.1,'--',HandleVisibility='off');
end
hold off
set(gca,'YDir','reverse');
grid minor;
legend(FontSize=15);
ylabel('P [dbar]',FontSize=15);

clr = 1:1:length(pCtd);
subplot(1,2,2)
plot(skLogn,kuLogn,'DisplayName','Lognormal','Color','#1f78b4',LineStyle='--',LineWidth=1.7);
hold on
plot(skLognN,kuLognN,'Color','#1f78b4',LineStyle='--',LineWidth=1.7,HandleVisibility='off');
for i = 1:n
    if strcmp(hypTest,"ks")
        if ks(2,i) < alphaKs
            plot(sk(i),ku(i),Marker="o",Color='k',HandleVisibility='off',MarkerSize=6);
        else
            plot(sk(i),ku(i),Marker="o",Color=[0.8 0.8 0.8],HandleVisibility='off');
        end
    else
        if ad(i) < alphaKs
            plot(sk(i),ku(i),Marker="o",Color='k',HandleVisibility='off',MarkerSize=6);
        else
            plot(sk(i),ku(i),Marker="o",Color=[0.8 0.8 0.8],HandleVisibility='off');
        end
    end
end
scatter(sk,ku,24,clr,"filled","o",HandleVisibility="off");
hold off
colormap(gca,cbrewer2("RdYlBu"));
cbar = colorbar;
cbar.Direction = "reverse";
cbar.Ticks = 1:10:length(pCtd);
cbar.TickLabels = pCtd(1):20:pCtd(101);
% cbar.Label.String = "P [dbar]";
ylim([1 kurtLimB]); xlim([skewLimA skewLimB]);
ylabel('Kurtosis',FontSize=15); xlabel('Skewness','FontSize',15);
legend(fontsize=15);
grid minor;

end

