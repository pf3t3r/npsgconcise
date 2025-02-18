function [] = plotKsSig(obsPerBin,threshold,limits,sigLimits,ks,trange,Sk,Ku)

subplot(1,3,1)
barh(obsPerBin,'FaceColor','#a6cee3');
hold on
xline(threshold);
hold off
set(gca,'XDir','reverse'); 
ylim(limits); set(gca,'YDir','reverse');
ylabel('$\Delta \sigma_0 [\textrm{kg m}^{-3}]$','Interpreter','latex');
set(gca,"YTick",5:5:length(trange),"YTickLabel",trange(5:5:length(trange)));
title('No. of observations');

subplot(1,3,2)
plot(ks(1,:),trange,'o-','Color','#a6cee3','DisplayName','Normal','LineWidth',1.5,'MarkerSize',5);
hold on
plot(ks(2,:),trange,'+--','Color','#1f78b4','DisplayName','Lognormal','LineWidth',1.5,'MarkerSize',5);
plot(ks(3,:),trange,'x-','Color','#b2df8a','DisplayName','Weibull','LineWidth',1.5,'MarkerSize',5);
plot(ks(4,:),trange,'.--','Color','#33a02c','DisplayName','Gamma','LineWidth',1.5,'MarkerSize',5);
hold off
grid minor;
xlabel('p-value');
ylim(sigLimits); set(gca,'YDir','reverse');
legend('Location','best');
yticklabels({});
title('KS Test');

subplot(1,3,3)
yyaxis left
plot(Sk,trange,'DisplayName','Skewness'); hold on
ylim(sigLimits); set(gca,'YDir','reverse'); yticklabels({});
yyaxis right
plot(Ku,trange,'DisplayName','Kurtosis');
ylim([-1.9 1.8]); set(gca,'YDir','reverse'); set(gca,'YTickLabel',{trange(5:5:length(trange))},'YColor','Black')
xline(3,'.','Mesokurtic','HandleVisibility','off');
xline(2.5,':','HandleVisibility','off');
xline(3.5,':','HandleVisibility','off');
xline(0,'.','Symmetric','HandleVisibility','off');
xline(-0.5,':','HandleVisibility','off');
xline(0.5,':','HandleVisibility','off');
hold off
grid minor;
legend('Location','south');
title('Moments');

end