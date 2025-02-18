function funcVisualiseDistributions(X,titleName)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

birS = fitdist(X,"BirnbaumSaunders");
expo = fitdist(X,"Exponential");
extV = fitdist(X,"ExtremeValue");
gamm = fitdist(X,"Gamma");
genE = fitdist(X,"GeneralizedExtremeValue");
genP = fitdist(X,"GeneralizedPareto");
halN = fitdist(X,"HalfNormal");
invG = fitdist(X,"InverseGaussian");
logi = fitdist(X,"Logistic");
logl = fitdist(X,"Loglogistic");
logn = fitdist(X,"Lognormal");
naka = fitdist(X,"Nakagami");
norm = fitdist(X,"Normal");
pois = fitdist(X,"Poisson");
rayl = fitdist(X,"Poisson");
rici = fitdist(X,"Rician");
tLoc = fitdist(X,"tLocationScale");
weib = fitdist(X,"Weibull");

x = 0:2*max(X);

yBirs = pdf(birS,x); yExpo = pdf(expo,x); yExtv = pdf(extV,x);
yGamm = pdf(gamm,x); yGene = pdf(genE,x); yGenp = pdf(genP,x);
yHaln = pdf(halN,x); yInvg = pdf(invG,x); yLogi = pdf(logi,x);
yLogl = pdf(logl,x); yLogn = pdf(logn,x); yNaka = pdf(naka,x);
yNorm = pdf(norm,x); yPois = pdf(pois,x); yRayl = pdf(rayl,x);
yRici = pdf(rici,x); yTloc = pdf(tLoc,x); yWeib = pdf(weib,x);

figure
yyaxis left
histogram(X);
hold on
yyaxis right
plot(x,yBirs,"DisplayName",'Birnbaum Saunders','Color','#a6cee3','LineStyle','-','Marker','none');
% plot(x,yExpo,"DisplayName",'Exponential');
plot(x,yExtv,"DisplayName",'Extreme Value','Color','#1f78b4','LineStyle','-','Marker','none');
plot(x,yGamm,"DisplayName",'Gamma','Color','#b2df8a','LineStyle','-','Marker','none','LineWidth',3);
plot(x,yGene,"DisplayName",'Generalized Extreme Value','Color','#33a02c','LineStyle','-','Marker','none');
% plot(x,yGenp,"DisplayName",'Generalized Pareto');
% plot(x,yHaln,"DisplayName",'Half Normal');
plot(x,yInvg,"DisplayName",'Inverse Gaussian','Color','#fb9a99','LineStyle','-','Marker','none');
plot(x,yLogi,"DisplayName",'Logistic','Color','#e31a1c','LineStyle','-','Marker','none');
plot(x,yLogl,"DisplayName",'LogLogistic','Color','#fdbf6f','LineStyle','-','Marker','none');
plot(x,yLogn,"DisplayName",'LogNormal','Color','#ff7f00','LineStyle','-','Marker','none','LineWidth',3);
plot(x,yNaka,"DisplayName",'Nakagami','Color','#cab2d6','LineStyle','-','Marker','none');
plot(x,yNorm,"DisplayName",'Normal','Color','#6a3d9a','LineStyle','-','Marker','none','LineWidth',3);
% plot(x,yPois,"DisplayName",'Poisson');
% plot(x,yRayl,"DisplayName",'Rayleigh');
plot(x,yRici,"DisplayName",'Rician','Color','#ffff99','LineStyle','-','Marker','none');
plot(x,yTloc,"DisplayName",'t Location-Scale','Color','#b15928','LineStyle','-','Marker','none');
plot(x,yWeib,"DisplayName",'Weibull','Color','#a6cee3','LineStyle','-','Marker','none','LineWidth',3);
hold off
legend();
title(titleName);

% figure;
% % histogram(chl_hplc);
% % hold on
% plot(birS);
% hold on
% % plot(expo); plot(extV); 
% plot(gamm); plot(genE);
% % plot(genP);
% % plot(halN);
% plot(invG); 
% % plot(logi); 
% plot(logl); plot(logn);
% % plot(naka);
% plot(norm);
% % plot(pois); plot(rayl);
% plot(rici); plot(tLoc); plot(weib);
% hold off
% legend();

end