function [a,b] = plotKs2(tr,ks,obs,sk,ku,limits,threshold,vuongRes,idObs,pV,hypTest,ad,testSel,logAxis)
%plotKs2
% INPUT: 
% OUTPUT: 
if nargin < 14
    logAxis = false;
end
if nargin < 13
    testSel = 4;
end
if nargin < 11
    hypTest = "ks";
end
if nargin < 8
    vuongRes = 0;
end
if nargin < 7
    threshold = 50;
end

% limits = [lim1 lim2];
% limits = limits;
tix = limits(1):10:limits(2);
disp(tix);
a = idObs(1); b = idObs(2);
n = length(tr);

% Create Annotations for Vuong's Test Results
annot = strings(1,n);
anClr = strings(1,n);
anClr(cellfun(@isempty,anClr)) = '#FFFFFF';
tmpEmph = strings(1,n); tmpEmph(cellfun(@isempty,tmpEmph)) = 'bold';

alphaHy = 0.005;
alphaLlr = 0.1;

vuongRes2 = nan(length(vuongRes),1);
if b > length(vuongRes)
    b = length(vuongRes);
end
vuongRes2(a:b) = vuongRes(a:b);
vuongRes = vuongRes2;

% 4.a. Vuong: Normal vs Lognormal vs Weibull vs Gamma
if testSel == 4
    for i = 1:n
        if strcmp(hypTest,"ks")
            if vuongRes(i) == 1 && ks(1,i) > alphaHy
                % Remove label if only one dist is not rejected by K-S.
                if length(find(ks(:,i)>alphaHy)) == 1
                    tmp = "";
                else
                    tmp = "Normal";
                end
                anClr(i) = '#a6cee3';
                if pV(1,i) > alphaLlr && ks(2,i) > alphaHy
                    tmpEmph(i) = 'normal';
                    tmp = append(tmp," L");
                end
                if pV(2,i) > alphaLlr && ks(3,i) > alphaHy
                    tmpEmph(i) = 'normal';
                    tmp = append(tmp," W");
                end
                if pV(3,i) > alphaLlr && ks(4,i) > alphaHy
                    tmpEmph(i) = 'normal';
                    tmp = append(tmp," G");
                end
                annot(i) = tmp;
            elseif vuongRes(i) == 2 && ks(2,i) > alphaHy
                % Remove label if only one dist is not rejected by K-S.
                if length(find(ks(:,i)>alphaHy)) == 1
                    tmp = "";
                else
                    tmp = "Lognormal";
                end
                anClr(i) = '#1f78b4';
                if pV(1,i) > alphaLlr && ks(1,i) > alphaHy
                    tmpEmph(i) = 'normal';
                    tmp = append(tmp," N");
                end
                if pV(5,i) > alphaLlr && ks(3,i) > alphaHy
                    tmpEmph(i) = 'normal';
                    tmp = append(tmp," W");
                end
                if pV(6,i) > alphaLlr && ks(4,i) > alphaHy
                    tmpEmph(i) = 'normal';
                    tmp = append(tmp," G");
                end
                annot(i) = tmp;
            elseif vuongRes(i) == 3 && ks(3,i) > alphaHy
                % Remove label if only one dist is not rejected by K-S.
                if length(find(ks(:,i)>alphaHy)) == 1
                    tmp = "";
                else
                    tmp = "Weibull";
                end
                anClr(i) = '#b2df8a';
                if pV(2,i) > alphaLlr && ks(1,i) > alphaHy
                    tmpEmph(i) = 'normal';
                    tmp = append(tmp," N");
                end
                if pV(5,i) > alphaLlr && ks(2,i) > alphaHy
                    tmpEmph(i) = 'normal';
                    tmp = append(tmp," L");
                end
                if pV(8,i) > alphaLlr && ks(4,i) > alphaHy
                    tmpEmph(i) = 'normal';
                    tmp = append(tmp," G");
                end
                annot(i) = tmp;
            elseif vuongRes(i) == 4 && ks(4,i) > alphaHy
                % Remove label if only one dist is not rejected by K-S.
                if length(find(ks(:,i)>alphaHy)) == 1
                    tmp = "";
                else
                    tmp = "Gamma";
                end
                anClr(i) = '#33a02c';
                if pV(6,i) > alphaLlr && ks(2,i) > alphaHy
                    tmpEmph(i) = 'normal';
                    tmp = append(tmp," L");
                end
                if pV(3,i) > alphaLlr && ks(1,i) > alphaHy
                    tmpEmph(i) = 'normal';
                    tmp = append(tmp," N");
                end
                if pV(8,i) > alphaLlr && ks(3,i) > alphaHy
                    tmpEmph(i) = 'normal';
                    tmp = append(tmp," W");
                end
                annot(i) = tmp;
            elseif vuongRes(i) == 0
                annot(i) = "";
            end
        elseif strcmp(hypTest,"ad")
            if vuongRes(i) == 1 && ad(1,i) > alphaHy
                % Remove label if only one dist is not rejected by K-S.
                if length(find(ad(:,i)>alphaHy)) == 1
                    tmp = "";
                else
                    tmp = "Normal";
                end
                anClr(i) = '#a6cee3';
                if pV(1,i) > alphaLlr && ad(2,i) > alphaHy
                    tmpEmph(i) = 'normal';
                    tmp = append(tmp," L");
                end
                if pV(2,i) > alphaLlr && ad(3,i) > alphaHy
                    tmpEmph(i) = 'normal';
                    tmp = append(tmp," W");
                end
                if pV(3,i) > alphaLlr && ad(4,i) > alphaHy
                    tmpEmph(i) = 'normal';
                    tmp = append(tmp," G");
                end
                annot(i) = tmp;
            elseif vuongRes(i) == 2 && ad(2,i) > alphaHy
                % Remove label if only one dist is not rejected by K-S.
                if length(find(ad(:,i)>alphaHy)) == 1
                    tmp = "";
                else
                    tmp = "Lognormal";
                end
                anClr(i) = '#1f78b4';
                if pV(1,i) > alphaLlr && ad(1,i) > alphaHy
                    tmpEmph(i) = 'normal';
                    tmp = append(tmp," N");
                end
                if pV(5,i) > alphaLlr && ad(3,i) > alphaHy
                    tmpEmph(i) = 'normal';
                    tmp = append(tmp," W");
                end
                if pV(6,i) > alphaLlr && ad(4,i) > alphaHy
                    tmpEmph(i) = 'normal';
                    tmp = append(tmp," G");
                end
                annot(i) = tmp;
            elseif vuongRes(i) == 3 && ad(3,i) > alphaHy
                % Remove label if only one dist is not rejected by K-S.
                if length(find(ad(:,i)>alphaHy)) == 1
                    tmp = "";
                else
                    tmp = "Weibull";
                end
                anClr(i) = '#b2df8a';
                if pV(2,i) > alphaLlr && ad(1,i) > alphaHy
                    tmpEmph(i) = 'normal';
                    tmp = append(tmp," N");
                end
                if pV(5,i) > alphaLlr && ad(2,i) > alphaHy
                    tmpEmph(i) = 'normal';
                    tmp = append(tmp," L");
                end
                if pV(8,i) > alphaLlr && ad(4,i) > alphaHy
                    tmpEmph(i) = 'normal';
                    tmp = append(tmp," G");
                end
                annot(i) = tmp;
            elseif vuongRes(i) == 4 && ad(4,i) > alphaHy
                % Remove label if only one dist is not rejected by K-S.
                if length(find(ad(:,i)>alphaHy)) == 1
                    tmp = "";
                else
                    tmp = "Gamma";
                end
                anClr(i) = '#33a02c';
                if pV(6,i) > alphaLlr && ad(2,i) > alphaHy
                    tmpEmph(i) = 'normal';
                    tmp = append(tmp," L");
                end
                if pV(3,i) > alphaLlr && ad(1,i) > alphaHy
                    tmpEmph(i) = 'normal';
                    tmp = append(tmp," N");
                end
                if pV(8,i) > alphaLlr && ad(3,i) > alphaHy
                    tmpEmph(i) = 'normal';
                    tmp = append(tmp," W");
                end
                annot(i) = tmp;
            elseif vuongRes(i) == 0
                annot(i) = "";
            end
        end
    end
elseif testSel == 2
    % 4.b. Vuong: Normal Vs. Lognormal Only
    if hypTest == "ks"
        for i = 1:n
            if vuongRes(i) == 1 && ks(1,i) > alphaHy
                annot(i) = "Normal";
                anClr(i) = '#c51b7d';
                if pV(1,i) > alphaLlr
                    tmpEmph(i) = 'normal';
                end
            elseif vuongRes(i) == 2 && ks(2,i) > alphaHy
                annot(i) = "Lognormal";
                anClr(i) = '#4d9221';
                if pV(1,i) > alphaLlr
                    tmpEmph(i) = 'normal';
                end
            else
                annot(i) = "";
            end
        end
    elseif hypTest == "ad"
        for i = 1:n
            if vuongRes(i) == 1 && ad(1,i) > alphaHy
                annot(i) = "Normal";
                anClr(i) = '#c51b7d';
                if pV(1,i) > alphaLlr
                    tmpEmph(i) = 'normal';
                end
            elseif vuongRes(i) == 2 && ad(2,i) > alphaHy
                annot(i) = "Lognormal";
                anClr(i) = '#4d9221';
                if pV(1,i) > alphaLlr
                    tmpEmph(i) = 'normal';
                end
            else
                annot(i) = "";
            end
        end
    end
end

% Lognormal family: generate theoretical skewness and kurtosis
sigTh = linspace(0,1,1000);
for i = 1:length(sigTh)
    skLogn(i) = (exp(sigTh(i)^2) + 2)*(sqrt(exp(sigTh(i)^2) - 1));
    kuLogn(i) = exp(4*sigTh(i)^2) + 2*exp(3*sigTh(i)^2) + 3*exp(2*sigTh(i)^2) - 3;
end

% % Gamma family: generate theoretical skewness and kurtosis
% % kTh = linspace(0.2,5000,10000);
% kTh = linspace(0.04,3000,1500000);
% for i = 1:length(kTh)
%     skGam(i) = 2/sqrt(kTh(i));
%     kuGam(i) = 6/kTh(i) + 3;
% end

% % Weibull family: generate theoretical skewness and kurtosis
% % kWbl = linspace(0,5,10000);
% kWbl = linspace(0.1,3.5,10000);
% for i = 1:length(kWbl)
%     skWbl(i) = ( gamma(1 + 3/kWbl(i)) - 3*gamma(1 + 1/kWbl(i))*gamma(1 + 2/kWbl(i)) + 2*(gamma(1 + 1/kWbl(i)))^3 ) ./ ...
%         ( gamma(1 + 2/kWbl(i)) -  (gamma(1 + 1/kWbl(i)))^2 )^(3/2);
%     kuWbl(i) = ( gamma(1 + 4/kWbl(i)) - 4*gamma(1 + 1/kWbl(i))*gamma(1 + 3/kWbl(i)) + 6*( (gamma(1 + 1/kWbl(i)) )^2)*gamma(1 + 2/kWbl(i)) - 3*( (gamma(1 + 1/kWbl(i)))^4 ) ) ./ ...
%        ( gamma(1 + 2/kWbl(i)) - ( gamma(1 + 1/kWbl(i)) )^2 )^2;
% end

% Negative Distributions
skLognN = -skLogn;
kuLognN = kuLogn;
% skGamN= -skGam;
% kuGamN = kuGam;
% skWblN = -skWbl;
% kuWblN = kuWbl;

subplot(1,3,3)
barh(obs(a:b),'FaceColor','#d3d3d3');
hold on
xline(threshold);
hold off
set(gca,'YDir','reverse');
yticklabels({});
xlabel('No. of Obs.',Interpreter='latex',FontSize=13);
% if b-a > 14 % b/c of bug with yticklabel
%     set(gca,"YTickLabel",limits(1):5:limits(2));
% else
% set(gca,"YTickLabel",limits(1):10:limits(2),"YTick",1:1+b-a);
% end
ylim([1 1+b-a]);

subplot(1,3,[1 2])
xline(alphaHy,'-','\color{black}\alpha=0.005',LineWidth=1.5,Color="#808080",HandleVisibility="off",LabelOrientation="horizontal",LabelHorizontalAlignment="center",FontSize=13);
hold on
if strcmp(hypTest,"ks")
    if testSel == 4
        plot(ks(1,:),tr,'o-','Color','#a6cee3','DisplayName','Normal','LineWidth',1.5,'MarkerSize',5);
        plot(ks(2,:),tr,'+--','Color','#1f78b4','DisplayName','Lognormal','LineWidth',1.5,'MarkerSize',5);
        plot(ks(3,:),tr,'x-','Color','#b2df8a','DisplayName','Weibull','LineWidth',1.5,'MarkerSize',5);
        plot(ks(4,:),tr,'.--','Color','#33a02c','DisplayName','Gamma','LineWidth',1.5,'MarkerSize',5);
    elseif testSel == 2
        plot(ks(1,:),tr,'o-','Color','#a6cee3','DisplayName','Normal','LineWidth',1.5,'MarkerSize',5,'HandleVisibility','off');
        plot(ks(2,:),tr,'+--','Color','#1f78b4','DisplayName','Lognormal','LineWidth',1.5,'MarkerSize',5,'HandleVisibility','off');
    end
    xlabel('K-S $p$-value',Interpreter='latex',FontSize=15);
else
    if testSel == 4
        plot(ad(1,:),tr,'o-','Color','#c51b7d','DisplayName','Normal','LineWidth',1.5,'MarkerSize',5);
        plot(ad(2,:),tr,'o-','Color','#4d9221','DisplayName','Lognormal','LineWidth',1.5,'MarkerSize',5);
        plot(ad(3,:),tr,'o-','Color','#d3d3d3','DisplayName','Weibull','LineWidth',1.5,'MarkerSize',5);
        plot(ad(4,:),tr,'o-','Color','#000000','DisplayName','Gamma','LineWidth',1.5,'MarkerSize',5);
    elseif testSel == 2
        for i = 1:n
            if vuongRes(i) == 1 && ad(1,i) > alphaHy & pV(1,i) > alphaLlr && ad(2,i) > alphaHy
                plot(ad(1,i),tr(i),'square','Color','#c51b7d','MarkerSize',15,HandleVisibility='off');
            elseif vuongRes(i) == 1 && ad(1,i) > alphaHy & pV(1,i) < alphaLlr && ad(2,i) > alphaHy
                plot(ad(1,i),tr(i),'square','Color','#c51b7d','MarkerSize',15,'LineWidth',4,HandleVisibility='off');
            elseif vuongRes(i) == 2 && ad(2,i) > alphaHy & pV(1,i) > alphaLlr && ad(1,i) > alphaHy
                plot(ad(2,i),tr(i),'square','Color','#4d9221','MarkerSize',15,HandleVisibility='off');
            elseif vuongRes(i) == 2 && ad(2,i) > alphaHy & pV(1,i) < alphaLlr && ad(1,i) > alphaHy
                plot(ad(2,i),tr(i),'square','Color','#4d9221','MarkerSize',15,'LineWidth',4,HandleVisibility='off');
            end
        end
        plot(ad(1,:),tr,'o-','Color','#c51b7d','DisplayName','Normal','LineWidth',1.5,'MarkerSize',5,'HandleVisibility','off');
        plot(ad(2,:),tr,'o-','Color','#4d9221','DisplayName','Lognormal','LineWidth',1.5,'MarkerSize',5,'HandleVisibility','off');
    end
    xlabel('A-D $p$-value',Interpreter='latex',FontSize=13);
end
if logAxis == true
    set(gca, 'XScale', 'log');
    %xline(0.05,':',HandleVisibility='off',LineWidth=1);
    %xline(0.1,':',HandleVisibility='off',LineWidth=1);
end
hold off
grid minor;
ylim([limits(1) limits(2)]); xlim([0.1*alphaHy 1]);
set(gca,'YDir','reverse');
set(gca,"YTick",limits(1):10:limits(2));
% yticklabels({});
ylabel('Pressure [dbar]',Interpreter='latex',FontSize=13);
% legend('Location','best',FontSize=11);
% title('K-S Test');

% zzs = 0.1*ones(n,1);
% subplot(1,6,4)
% for i = 1:n
%     text(zzs(i),tr(i),annot(i),FontSize=10,Color=anClr(i),FontWeight=tmpEmph(i));
% end
% % % For Normal-Lognormal Comparison ONLY
% % hold on
% % pV(1,obs<threshold) = nan;
% % for i = 1:n
% %     if pV(1,i) > 0.01
% %         scatter(pV(1,i),tr(i),[],"black");
% %     end
% % end
% % hold off
% grid minor;
% ylim([limits(1) limits(2)]); 
% set(gca,'YDir','reverse');
% set(gca,"YTick",limits(1):10:limits(2));
% yticklabels({});
% xticklabels({' ' ,' '});
% % set(gca,"XTick",[]);
% xlabel('Vuong LLR',FontSize=15);


% subplot(1,4,3)
% yyaxis left
% plot(sk,tr,'DisplayName','Skewness'); hold on
% ylim([limits(1) limits(2)]);
% set(gca,'YDir','reverse');
% % yticklabels({});
% yyaxis right
% plot(ku,tr,'DisplayName','Kurtosis');
% ylim([limits(1) limits(2)]); 
% set(gca,'YDir','reverse');
% % set(gca,'YTickLabel',{tr(5:5:n)},'YColor','Black')
% xline(3,'.','Mesokurtic','HandleVisibility','off');
% xline(2.5,':','HandleVisibility','off');
% xline(3.5,':','HandleVisibility','off');
% xline(0,'.','Symmetric','HandleVisibility','off');
% xline(-0.5,':','HandleVisibility','off');
% xline(0.5,':','HandleVisibility','off');
% hold off
% grid minor;
% legend('Location','south');
% title('Moments');

% % SK-KU plot: output only skewness/kurtosis values at depths defined by
% % y-limits in function.
% tmp = [];
% for i = a:b
%     if ~isnan(sum(ks(:,i)))
%         tmp = [tmp i];
%     end
% end
% tr2 = tr(tmp);
% sk2 = sk(tmp);
% ku2 = ku(tmp);
% clear tmp;
% 
% kurtLimB = 10; skewLimA = 0; skewLimB = 2.5;
% if max(ku2) > 10 & min(sk2) < 0
%     kurtLimB = max(ku2) + 1;
%     skewLimA = min(sk2) - 0.1;
%     skewLimB = max(sk2) + 0.1;
% elseif max(ku2) > 10
%     kurtLimB = max(ku2) + 1;
%     skewLimB = max(sk2) + 0.1;
% elseif min(sk2) < 0 
%     skewLimA = min(sk2) - 0.1;
% elseif max(sk2) > 2.5
%     kurtLimB = max(ku2) + 1;
%     skewLimB = max(sk2) + 0.1;
% else 
%     kurtLimB = 10;
%     skewLimA = 0;
%     skewLimB = 2.5;
% end
% 
% % % error bars for sk-ku
% % obsTmp = obs(obs>=50);
% % n2 = length(obsTmp);
% % yneg = nan(n2,1); ypos = nan(n2,1); xneg = nan(n2,1); xpos = nan(n2,1);
% % for i = 1:n2
% %     if obsTmp(i) > 300
% %         yneg(i) = -1.02; ypos(i) = 0.92;
% %         xneg(i) = -0.23; xpos(i) = 0.22;
% %     elseif obsTmp(i) > 250
% %         yneg(i) = -1.06; ypos(i) = 0.96;
% %         xneg(i) = -0.24; xpos(i) = 0.25;
% %     elseif obsTmp(i) > 200
% %         yneg(i) = -1.13; ypos(i) = 0.96;
% %         xneg(i) = -0.27; xpos(i) = 0.26;
% %     elseif obsTmp(i) > 150
% %         yneg(i) = -1.22; ypos(i) = 1.08;
% %         xneg(i) = -0.30; xpos(i) = 0.29;
% %     elseif obsTmp(i) > 100
% %         yneg(i) = -1.30; ypos(i) = 1.15;
% %         xneg(i) = -0.34; xpos(i) = 0.33;
% %     else
% %         yneg(i) = -1.24; ypos(i) = 1.11;
% %         xneg(i) = -0.51; xpos(i) = 0.38;
% %     end
% % end
% 
% subplot(1,6,[5.1 6])
% scatter(0,3,[],[0.6509803921568628 0.807843137254902 0.8901960784313725],'DisplayName','Norm.',Marker='o',LineWidth=3);
% hold on
% plot(skLogn,kuLogn,'DisplayName','Logn.','Color','#1f78b4',LineStyle='--',LineWidth=1.7);
% plot(skLognN,kuLognN,'Color','#1f78b4',LineStyle='--',LineWidth=1.7,HandleVisibility='off');
% plot(skWbl,kuWbl,'DisplayName','Weib.','Color','#b2df8a',LineStyle='-',LineWidth=1.7);
% plot(skWblN,kuWblN,'Color','#b2df8a',LineStyle='-',LineWidth=1.7,HandleVisibility='off');
% plot(skGam,kuGam,'DisplayName','Gam.','Color','#33a02c',LineStyle='--',LineWidth=1.7);
% plot(skGamN,kuGamN,'Color','#33a02c',LineStyle='--',LineWidth=1.7,HandleVisibility='off');
% scatter(2,9,'DisplayName','Exp.',Marker='+',LineWidth=1);
% scatter(0,9/5,'DisplayName','Uni.',Marker='o',LineWidth=1);
% scatter(0,21/5,'DisplayName','Logi.',Marker='.',LineWidth=1);
% scatter(1.1395,5.4,'DisplayName','LEV',Marker='x',LineWidth=1);
% % errorbar(sk2,ku2,yneg,ypos,xneg,xpos,'o','Color',[0.6 0.6 0.6],'HandleVisibility','off');
% scatter(sk2,ku2,[],[0.8 0.8 0.8],HandleVisibility='off');
% clr = 1:1:length(tr2);
% scatter(sk2,ku2,24,clr,"filled","o",HandleVisibility="off");
% % colormap(gca,flipud(colormap("hot")));
% colormap(gca,cbrewer2("RdYlBu"));
% cbar = colorbar;
% cbar.Direction = "reverse";
% cbar.Ticks = 1:1:length(tr2);
% % cbar.TickLabels = tr2(1):10:tr2(end);
% cbar.TickLabels = tr2;
% cbar.Label.String = "P [dbar]";
% hold off
% grid minor;
% ylim([1 kurtLimB]); xlim([skewLimA skewLimB]);
% xlabel('Skewness',FontSize=15); ylabel('Kurtosis',FontSize=15);
% lgd = legend('Location','best');
% title(lgd,'Distributions');
% % title('Skewness vs. Kurtosis');

end