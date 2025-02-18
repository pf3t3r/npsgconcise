% Skewness-kurtosis analysis for L2 bottle data.

% This script duplicates a lot of L2_bot right now. Must clean this up.

clear; clc; close all;
addpath("baroneRoutines\"); addpath("func\");
set(groot, "defaultFigureUnits", "centimeters", "defaultFigurePosition", [3 3 10 15]);

% Load maximum mixed-layer depth 'MLD' and cruise-averaged deep chlorophyll
% maximum depth 'DCM'.
mld = load("mldVals.mat").maxMld; % single maximum per cruise
dcm = load("output/dcm.mat").dcm; % pDcm + sigmaDcm (all casts, crn 1-329)

testSel = 2; % 2 = norm + logn; 4 = norm + logn + weib + gamm

tmp = importdata("data/L2/hplcChla_88-21_200.txt");

%%

id = num2str(tmp.data(:,1));
p = tmp.data(:,4);
c = tmp.data(:,5);
threshold = 30;

idOld = id;
% %%% seasonal analysis
% if season ~= 0
%     botidIn = tmp.data(:,2);
%     n = length(p);
%     botidIn(botidIn==-9) = nan;
%     botId2 = num2str(botidIn);
%     botMth = nan(n,1);
%     for i = 1:n
%         tmpX = str2num(botId2(i,1:end-4));
%         if ~isnan(tmpX)
%             botMth(i) = tmpX;
%         end
%     end
%     winter = nan(n,1); spring = nan(n,1); summer = nan(n,1); autumn = nan(n,1);
%     for i = 1:n
%         tmpY = botMth(i);
%         if (tmpY == 12) || (tmpY == 1) || (tmpY == 2)
%             winter(i) = 1;
%         end
%         if (tmpY == 3) || (tmpY == 4) || (tmpY == 5)
%             spring(i) = 1;
%         end
%         if (tmpY == 6) || (tmpY == 7) || (tmpY == 8)
%             summer(i) = 1;
%         end
%         if (tmpY == 9) || (tmpY == 10) || (tmpY == 11)
%             autumn(i) = 1;
%         end
%     end
% 
%     winIds = []; sprIds = []; sumIds = []; autIds = []; 
%     for i = 1:n
%         if winter(i) == 1
%             winIds = [winIds i];
%         end
%         if spring(i) == 1
%             sprIds = [sprIds i];
%         end
%         if summer(i) == 1
%             sumIds = [sumIds i];
%         end
%         if autumn(i) == 1
%             autIds = [autIds i];
%         end
%     end
% 
%     if season == 1
%         c = c(winIds);
%         p = p(winIds);
%         id = id(winIds,:);
%         %nSea = length(winIds);
%     elseif season == 2
%         c = c(sprIds);
%         p = p(sprIds);
%         id = id(sprIds,:);
%         %nSea = length(sprIds);
%     elseif season == 3
%         c = c(sumIds);
%         p = p(sumIds);
%         id = id(sumIds,:);
%         %nSea = length(sumIds);
%     elseif season == 4
%         c = c(autIds);
%         p = p(autIds);
%         id = id(autIds,:);
%         %nSea = length(autIds);
%     end
% end
% %%% end seasonal analysis

% 2. Extract data beneath ML
[idSubml,pSubml,cSubml] = extractSMLC(id,p,c,mld);

% 3. Calculate KS p-value, skewness, kurtosis
% ..., centre around DCM (?)
[pr,ks,obs,sk,ku,rV,pV,ad,X_out,bA] = ksOfLagrangian(idSubml,pSubml,dcm,cSubml,threshold);

% 3.a. Intercomparison of results from Vuong's Test: easily see best
% distribution at each depth.
vuongRes = zeros(1,length(pr));
rV(isnan(rV)) = 0;

if testSel==4
    for i = 1:length(pr)
        %disp(i);
        if rV(1,i) & rV(2,i) & rV(3,i) > 0
            %disp('Normal');
            vuongRes(i) = 1;
        elseif rV(1,i) < 0 & rV(5,i) > 0 & rV(6,i) > 0
            %disp('Lognormal');
            vuongRes(i) = 2;
        elseif rV(2,i) < 0 & rV(5,i) < 0 & rV(8,i) > 0
            %disp('Weibull');
            vuongRes(i) = 3;
        elseif rV(3,i) < 0 & rV(6,i) < 0 & rV(8,i) < 0
            %disp('Gamma');
            vuongRes(i) = 4;
        end
    end
    rV(rV==0) = nan;
elseif testSel == 2
    for i = 1:length(pr)
        if rV(1,i)  > 0
            vuongRes(i) = 1;
        elseif rV(1,i) < 0
            vuongRes(i) = 2;
        end
    end
    rV(rV==0) = nan;
end

% limits = [pr(barchartLimits(1)) pr(barchartLimits(2))];
limits = [-60 60];
obsId = [7 19];

% 4. Plot results
ax = figure;

% % limits = [lim1 lim2];
% % limits = limits;
% tix = limits(1):10:limits(2);
% disp(tix);
a = obsId(1); b = obsId(2);
n = length(pr);

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

% % 4.a. Vuong: Normal vs Lognormal vs Weibull vs Gamma
% if testSel == 4
%     for i = 1:n
%         if strcmp(hypTest,"ks")
%             if vuongRes(i) == 1 && ks(1,i) > alphaHy
%                 % Remove label if only one dist is not rejected by K-S.
%                 if length(find(ks(:,i)>alphaHy)) == 1
%                     tmp = "";
%                 else
%                     tmp = "Normal";
%                 end
%                 anClr(i) = '#a6cee3';
%                 if pV(1,i) > alphaLlr && ks(2,i) > alphaHy
%                     tmpEmph(i) = 'normal';
%                     tmp = append(tmp," L");
%                 end
%                 if pV(2,i) > alphaLlr && ks(3,i) > alphaHy
%                     tmpEmph(i) = 'normal';
%                     tmp = append(tmp," W");
%                 end
%                 if pV(3,i) > alphaLlr && ks(4,i) > alphaHy
%                     tmpEmph(i) = 'normal';
%                     tmp = append(tmp," G");
%                 end
%                 annot(i) = tmp;
%             elseif vuongRes(i) == 2 && ks(2,i) > alphaHy
%                 % Remove label if only one dist is not rejected by K-S.
%                 if length(find(ks(:,i)>alphaHy)) == 1
%                     tmp = "";
%                 else
%                     tmp = "Lognormal";
%                 end
%                 anClr(i) = '#1f78b4';
%                 if pV(1,i) > alphaLlr && ks(1,i) > alphaHy
%                     tmpEmph(i) = 'normal';
%                     tmp = append(tmp," N");
%                 end
%                 if pV(5,i) > alphaLlr && ks(3,i) > alphaHy
%                     tmpEmph(i) = 'normal';
%                     tmp = append(tmp," W");
%                 end
%                 if pV(6,i) > alphaLlr && ks(4,i) > alphaHy
%                     tmpEmph(i) = 'normal';
%                     tmp = append(tmp," G");
%                 end
%                 annot(i) = tmp;
%             elseif vuongRes(i) == 3 && ks(3,i) > alphaHy
%                 % Remove label if only one dist is not rejected by K-S.
%                 if length(find(ks(:,i)>alphaHy)) == 1
%                     tmp = "";
%                 else
%                     tmp = "Weibull";
%                 end
%                 anClr(i) = '#b2df8a';
%                 if pV(2,i) > alphaLlr && ks(1,i) > alphaHy
%                     tmpEmph(i) = 'normal';
%                     tmp = append(tmp," N");
%                 end
%                 if pV(5,i) > alphaLlr && ks(2,i) > alphaHy
%                     tmpEmph(i) = 'normal';
%                     tmp = append(tmp," L");
%                 end
%                 if pV(8,i) > alphaLlr && ks(4,i) > alphaHy
%                     tmpEmph(i) = 'normal';
%                     tmp = append(tmp," G");
%                 end
%                 annot(i) = tmp;
%             elseif vuongRes(i) == 4 && ks(4,i) > alphaHy
%                 % Remove label if only one dist is not rejected by K-S.
%                 if length(find(ks(:,i)>alphaHy)) == 1
%                     tmp = "";
%                 else
%                     tmp = "Gamma";
%                 end
%                 anClr(i) = '#33a02c';
%                 if pV(6,i) > alphaLlr && ks(2,i) > alphaHy
%                     tmpEmph(i) = 'normal';
%                     tmp = append(tmp," L");
%                 end
%                 if pV(3,i) > alphaLlr && ks(1,i) > alphaHy
%                     tmpEmph(i) = 'normal';
%                     tmp = append(tmp," N");
%                 end
%                 if pV(8,i) > alphaLlr && ks(3,i) > alphaHy
%                     tmpEmph(i) = 'normal';
%                     tmp = append(tmp," W");
%                 end
%                 annot(i) = tmp;
%             elseif vuongRes(i) == 0
%                 annot(i) = "";
%             end
%         elseif strcmp(hypTest,"ad")
%             if vuongRes(i) == 1 && ad(1,i) > alphaHy
%                 % Remove label if only one dist is not rejected by K-S.
%                 if length(find(ad(:,i)>alphaHy)) == 1
%                     tmp = "";
%                 else
%                     tmp = "Normal";
%                 end
%                 anClr(i) = '#a6cee3';
%                 if pV(1,i) > alphaLlr && ad(2,i) > alphaHy
%                     tmpEmph(i) = 'normal';
%                     tmp = append(tmp," L");
%                 end
%                 if pV(2,i) > alphaLlr && ad(3,i) > alphaHy
%                     tmpEmph(i) = 'normal';
%                     tmp = append(tmp," W");
%                 end
%                 if pV(3,i) > alphaLlr && ad(4,i) > alphaHy
%                     tmpEmph(i) = 'normal';
%                     tmp = append(tmp," G");
%                 end
%                 annot(i) = tmp;
%             elseif vuongRes(i) == 2 && ad(2,i) > alphaHy
%                 % Remove label if only one dist is not rejected by K-S.
%                 if length(find(ad(:,i)>alphaHy)) == 1
%                     tmp = "";
%                 else
%                     tmp = "Lognormal";
%                 end
%                 anClr(i) = '#1f78b4';
%                 if pV(1,i) > alphaLlr && ad(1,i) > alphaHy
%                     tmpEmph(i) = 'normal';
%                     tmp = append(tmp," N");
%                 end
%                 if pV(5,i) > alphaLlr && ad(3,i) > alphaHy
%                     tmpEmph(i) = 'normal';
%                     tmp = append(tmp," W");
%                 end
%                 if pV(6,i) > alphaLlr && ad(4,i) > alphaHy
%                     tmpEmph(i) = 'normal';
%                     tmp = append(tmp," G");
%                 end
%                 annot(i) = tmp;
%             elseif vuongRes(i) == 3 && ad(3,i) > alphaHy
%                 % Remove label if only one dist is not rejected by K-S.
%                 if length(find(ad(:,i)>alphaHy)) == 1
%                     tmp = "";
%                 else
%                     tmp = "Weibull";
%                 end
%                 anClr(i) = '#b2df8a';
%                 if pV(2,i) > alphaLlr && ad(1,i) > alphaHy
%                     tmpEmph(i) = 'normal';
%                     tmp = append(tmp," N");
%                 end
%                 if pV(5,i) > alphaLlr && ad(2,i) > alphaHy
%                     tmpEmph(i) = 'normal';
%                     tmp = append(tmp," L");
%                 end
%                 if pV(8,i) > alphaLlr && ad(4,i) > alphaHy
%                     tmpEmph(i) = 'normal';
%                     tmp = append(tmp," G");
%                 end
%                 annot(i) = tmp;
%             elseif vuongRes(i) == 4 && ad(4,i) > alphaHy
%                 % Remove label if only one dist is not rejected by K-S.
%                 if length(find(ad(:,i)>alphaHy)) == 1
%                     tmp = "";
%                 else
%                     tmp = "Gamma";
%                 end
%                 anClr(i) = '#33a02c';
%                 if pV(6,i) > alphaLlr && ad(2,i) > alphaHy
%                     tmpEmph(i) = 'normal';
%                     tmp = append(tmp," L");
%                 end
%                 if pV(3,i) > alphaLlr && ad(1,i) > alphaHy
%                     tmpEmph(i) = 'normal';
%                     tmp = append(tmp," N");
%                 end
%                 if pV(8,i) > alphaLlr && ad(3,i) > alphaHy
%                     tmpEmph(i) = 'normal';
%                     tmp = append(tmp," W");
%                 end
%                 annot(i) = tmp;
%             elseif vuongRes(i) == 0
%                 annot(i) = "";
%             end
%         end
%     end
% elseif testSel == 2
%     % 4.b. Vuong: Normal Vs. Lognormal Only
%     if hypTest == "ks"
%         for i = 1:n
%             if vuongRes(i) == 1 && ks(1,i) > alphaHy
%                 annot(i) = "Normal";
%                 anClr(i) = '#c51b7d';
%                 if pV(1,i) > alphaLlr
%                     tmpEmph(i) = 'normal';
%                 end
%             elseif vuongRes(i) == 2 && ks(2,i) > alphaHy
%                 annot(i) = "Lognormal";
%                 anClr(i) = '#4d9221';
%                 if pV(1,i) > alphaLlr
%                     tmpEmph(i) = 'normal';
%                 end
%             else
%                 annot(i) = "";
%             end
%         end
%     elseif hypTest == "ad"
%         for i = 1:n
%             if vuongRes(i) == 1 && ad(1,i) > alphaHy
%                 annot(i) = "Normal";
%                 anClr(i) = '#c51b7d';
%                 if pV(1,i) > alphaLlr
%                     tmpEmph(i) = 'normal';
%                 end
%             elseif vuongRes(i) == 2 && ad(2,i) > alphaHy
%                 annot(i) = "Lognormal";
%                 anClr(i) = '#4d9221';
%                 if pV(1,i) > alphaLlr
%                     tmpEmph(i) = 'normal';
%                 end
%             else
%                 annot(i) = "";
%             end
%         end
%     end
% end

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

% subplot(1,3,3)
% barh(obs(a:b),'FaceColor','#d3d3d3');
% hold on
% xline(threshold);
% hold off
% set(gca,'YDir','reverse');
% % ylabel('Pressure [dbar]',FontSize=15);
% yticklabels({});
% xlabel('No. of Observations',Interpreter='latex',FontSize=13);
% % if b-a > 14 % b/c of bug with yticklabel
% %     set(gca,"YTickLabel",limits(1):5:limits(2));
% % else
% set(gca,"YTickLabel",limits(1):10:limits(2),"YTick",1:1+b-a);
% % end
% ylim([1 1+b-a]);

% subplot(1,3,[1 2])
% xline(alphaHy,DisplayName='\alpha',LineWidth=1.5,Color="#808080");
% hold on
% if strcmp(hypTest,"ks")
%     if testSel == 4
%         plot(ks(1,:),pr,'o-','Color','#a6cee3','DisplayName','Normal','LineWidth',1.5,'MarkerSize',5);
%         plot(ks(2,:),pr,'+--','Color','#1f78b4','DisplayName','Lognormal','LineWidth',1.5,'MarkerSize',5);
%         plot(ks(3,:),pr,'x-','Color','#b2df8a','DisplayName','Weibull','LineWidth',1.5,'MarkerSize',5);
%         plot(ks(4,:),pr,'.--','Color','#33a02c','DisplayName','Gamma','LineWidth',1.5,'MarkerSize',5);
%     elseif testSel == 2
%         plot(ks(1,:),pr,'o-','Color','#a6cee3','DisplayName','Normal','LineWidth',1.5,'MarkerSize',5);
%         plot(ks(2,:),pr,'+--','Color','#1f78b4','DisplayName','Lognormal','LineWidth',1.5,'MarkerSize',5);
%     end
%     xlabel('K-S $p$-value',Interpreter='latex',FontSize=15);
% else
%     if testSel == 4
%         plot(ad(1,:),pr,'o-','Color','#a6cee3','DisplayName','Normal','LineWidth',1.5,'MarkerSize',5);
%         plot(ad(2,:),pr,'+--','Color','#1f78b4','DisplayName','Lognormal','LineWidth',1.5,'MarkerSize',5);
%         plot(ad(3,:),pr,'x-','Color','#b2df8a','DisplayName','Weibull','LineWidth',1.5,'MarkerSize',5);
%         plot(ad(4,:),pr,'.--','Color','#33a02c','DisplayName','Gamma','LineWidth',1.5,'MarkerSize',5);
%     elseif testSel == 2
%         for i = 1:n
%             if vuongRes(i) == 1 && ad(1,i) > alphaHy & pV(1,i) > alphaLlr
%                 plot(ad(1,i),pr(i),'square','Color','#c51b7d','MarkerSize',15,HandleVisibility='off');
%             elseif vuongRes(i) == 1 && ad(1,i) > alphaHy & pV(1,i) < alphaLlr
%                 plot(ad(1,i),pr(i),'square','Color','#c51b7d','MarkerSize',15,'LineWidth',4,HandleVisibility='off');
%             elseif vuongRes(i) == 2 && ad(2,i) > alphaHy & pV(1,i) > alphaLlr
%                 plot(ad(2,i),pr(i),'square','Color','#4d9221','MarkerSize',15,HandleVisibility='off');
%             elseif vuongRes(i) == 2 && ad(2,i) > alphaHy & pV(1,i) < alphaLlr
%                 plot(ad(2,i),pr(i),'square','Color','#4d9221','MarkerSize',15,'LineWidth',4,HandleVisibility='off');
%             end
%         end
%         plot(ad(1,:),pr,'o-','Color','#c51b7d','DisplayName','Normal','LineWidth',1.5,'MarkerSize',5);
%         plot(ad(2,:),pr,'+--','Color','#4d9221','DisplayName','Lognormal','LineWidth',1.5,'MarkerSize',5);
%     end
%     xlabel('A-D $p$-value + Vuong Ratio',Interpreter='latex',FontSize=13);
% end
% if logAxis == true
%     set(gca, 'XScale', 'log');
%     xline(0.05,':',HandleVisibility='off',LineWidth=1);
%     xline(0.1,':',HandleVisibility='off',LineWidth=1);
% end
% hold off
% grid minor;
% ylim([limits(1) limits(2)]); xlim([0.6*alphaHy 1]);
% set(gca,'YDir','reverse');
% set(gca,"YTick",limits(1):10:limits(2));
% % yticklabels({});
% ylabel('Pressure [dbar]',Interpreter='latex',FontSize=13);
% legend('Location','best',FontSize=11);
% % title('K-S Test');

% zzs = 0.1*ones(n,1);
% subplot(1,6,4)
% for i = 1:n
%     text(zzs(i),pr(i),annot(i),FontSize=10,Color=anClr(i),FontWeight=tmpEmph(i));
% end
% % % For Normal-Lognormal Comparison ONLY
% % hold on
% % pV(1,obs<threshold) = nan;
% % for i = 1:n
% %     if pV(1,i) > 0.01
% %         scatter(pV(1,i),pr(i),[],"black");
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
% plot(sk,pr,'DisplayName','Skewness'); hold on
% ylim([limits(1) limits(2)]);
% set(gca,'YDir','reverse');
% % yticklabels({});
% yyaxis right
% plot(ku,pr,'DisplayName','Kurtosis');
% ylim([limits(1) limits(2)]); 
% set(gca,'YDir','reverse');
% % set(gca,'YTickLabel',{pr(5:5:n)},'YColor','Black')
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

% SK-KU plot: output only skewness/kurtosis values at depths defined by
% y-limits in function.
tmp = [];
for i = a:b
    if ~isnan(sum(ks(:,i)))
        tmp = [tmp i];
    end
end
tr2 = pr(tmp);
sk2 = sk(tmp);
ku2 = ku(tmp);
clear tmp;

kurtLimB = 10; skewLimA = 0; skewLimB = 2.5;
if max(ku2) > 10 & min(sk2) < 0
    kurtLimB = max(ku2) + 1;
    skewLimA = min(sk2) - 0.1;
    skewLimB = max(sk2) + 0.1;
elseif max(ku2) > 10
    kurtLimB = max(ku2) + 1;
    skewLimB = max(sk2) + 0.1;
elseif min(sk2) < 0 
    skewLimA = min(sk2) - 0.1;
elseif max(sk2) > 2.5
    kurtLimB = max(ku2) + 1;
    skewLimB = max(sk2) + 0.1;
else 
    kurtLimB = 10;
    skewLimA = 0;
    skewLimB = 2.5;
end
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

scatter(0,3,72,[0.2 0.2 0.2],'DisplayName','Normal',Marker='pentagram',LineWidth=2.5);
hold on
plot(skLogn,kuLogn,'DisplayName','Lognormal','Color','#808080',LineStyle='-',LineWidth=1.3);
plot(skLognN,kuLognN,'Color','#808080',LineStyle='-',LineWidth=1.3,HandleVisibility='off');
% plot(skWbl,kuWbl,'DisplayName','Weib.','Color','#b2df8a',LineStyle='-',LineWidth=1.7);
% plot(skWblN,kuWblN,'Color','#b2df8a',LineStyle='-',LineWidth=1.7,HandleVisibility='off');
% plot(skGam,kuGam,'DisplayName','Gam.','Color','#33a02c',LineStyle='--',LineWidth=1.7);
% plot(skGamN,kuGamN,'Color','#33a02c',LineStyle='--',LineWidth=1.7,HandleVisibility='off');
% scatter(2,9,'DisplayName','Exp.',Marker='+',LineWidth=1);
% scatter(0,9/5,'DisplayName','Uni.',Marker='o',LineWidth=1);
% scatter(0,21/5,'DisplayName','Logi.',Marker='.',LineWidth=1);
% scatter(1.1395,5.4,'DisplayName','LEV',Marker='x',LineWidth=1);
% errorbar(sk2,ku2,yneg,ypos,xneg,xpos,'o','Color',[0.6 0.6 0.6],'HandleVisibility','off');
scatter(sk2,ku2,72,[0.8 0.8 0.8],HandleVisibility='off');
clr = 1:1:length(tr2);
scatter(sk2,ku2,54,clr,"filled","o",HandleVisibility="off");
% colormap(gca,flipud(colormap("hot")));
colormap(gca,flipud(cbrewer2("PiYG")));
cbar = colorbar;
cbar.Direction = "reverse";
cbar.Ticks = 1:1:length(tr2);
% cbar.TickLabels = tr2(1):10:tr2(end);
cbar.TickLabels = tr2;
cbar.Label.String = "P [dbar]";
cbar.Label.Position = [0.7 1-0.7];
cbar.Label.Rotation = 0;
% add polynomial
[skS,id] = sort(sk2);
kuS = ku2(id);
[p,S] = polyfit(skS,kuS,2);
[f,delta] = polyval(p,skS,S);
plot(skS,f,'r-',DisplayName="Fit");
plot(skS,f+2*delta,'m--',skS,f-2*delta,'m--');
hold off
grid minor;
ylim([1 kurtLimB]); xlim([skewLimA skewLimB]);
xlabel('Skewness',FontSize=13,Interpreter='latex'); 
% ylabel('Kurtosis',FontSize=13,Interpreter='latex');
yticklabels({});
% lgd = legend('Location','best');
% title(lgd,'Distributions');
% title('Skewness vs. Kurtosis');

%%
title("L2","Interpreter","latex",FontSize=13);
exportgraphics(ax,"figures/L2/bottle/log/chla_ad_skKu.png");