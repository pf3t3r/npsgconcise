% Skewness-kurtosis analysis for L1 bottle data.

% This script duplicates a lot of L1_bot right now. Must clean this up.

clear; clc; close all;
addpath("baroneRoutines\"); addpath("func\");
set(groot, "defaultFigureUnits", "centimeters", "defaultFigurePosition", [3 3 10 15]);

ctdData = importdata("datafiles\ctd_iso_ALL.mat").ctd;
maxMld = nan(329,1);
for i = 1:329
    if ~isnan([ctdData(i).mld003])
        maxMld(i) = max([ctdData(i).mld003]);
    end
end

clear ctdData i;
save mldVals.mat maxMld;

tmp = importdata("data/L1/hplcChla_88-21_150.txt");
pIn = tmp.data(:,4);
cIn = tmp.data(:,5);
idIn = tmp.data(:,1);

% %%% seasonal analysis
% if season ~= 0
%     botidIn = tmp.data(:,2);
%     n = length(pIn);
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
%         cIn = cIn(winIds);
%         pIn = pIn(winIds);
%         idIn = idIn(winIds);
%         %nSea = length(winIds);
%     elseif season == 2
%         cIn = cIn(sprIds);
%         pIn = pIn(sprIds);
%         idIn = idIn(sprIds);
%         %nSea = length(sprIds);
%     elseif season == 3
%         cIn = cIn(sumIds);
%         pIn = pIn(sumIds);
%         idIn = idIn(sumIds);
%         %nSea = length(sumIds);
%     elseif season == 4
%         cIn = cIn(autIds);
%         pIn = pIn(autIds);
%         idIn = idIn(autIds);
%         %nSea = length(autIds);
%     end
% end
% %%% end seasonal analysis

% 2. Extract data in ML
[idOut,pOut,cOut] = extractMldVals(idIn,pIn,cIn,maxMld);

% 3. Bin data
[~,pOutB,cOutB,~,~] = cleanAndBin(pOut,cOut,idOut');

% 4. Calculate KS p-value, skewness, kurtosis, Vuong Parameters
[ks,obs,p,sk,ku,rV,pV,ad] = ksOfBinnedCon(cOutB,pOutB,10,30);

% 4.a. Intercomparison of results from Vuong's Test: easily see best
% distribution at each depth.
vuongRes = nan(1,length(p));

testSel = 2;
% 4.a.i Default Case.
if testSel == 4
    for i = 1:length(p)
        if rV(1,i) & rV(2,i) & rV(3,i) > 0
            disp('Normal');
            vuongRes(i) = 1;
        elseif rV(1,i) < 0 & rV(5,i) > 0 & rV(6,i) > 0
            disp('Lognormal');
            vuongRes(i) = 2;
        elseif rV(2,i) < 0 & rV(5,i) < 0 & rV(8,i) > 0
            disp('Weibull');
            vuongRes(i) = 3;
        elseif rV(3,i) < 0 & rV(6,i) < 0 & rV(8,i) < 0
            disp('Gamma');
            vuongRes(i) = 4;
        end
    end
elseif testSel == 2
% 4.a.ii. Normal-Lognormal Case ONLY.
    for i = 1:length(p)
        if rV(1,i) > 0
            disp('Normal');
            vuongRes(i) = 1;
        elseif rV(1,i) < 0
            disp('Lognormal');
            vuongRes(i) = 2;
        end
    end
end

% 5. Plot results
ax = figure;

n = length(p);

% % Create Annotations for Vuong's Test Results
% annot = strings(1,n);
% anClr = strings(1,n);
% anClr(cellfun(@isempty,anClr)) = '#FFFFFF';
% tmpEmph = strings(1,n); tmpEmph(cellfun(@isempty,tmpEmph)) = 'bold';
% % Default Case
% 
% alphaHy = 0.005;
% alphaLlr = 0.10;

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
%             % A-D
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
% 
%         end
%     end
% elseif testSel == 2
%     % Normal-Lognormal Case
%     for i = 1:n
%         if strcmp(hypTest,"ks")
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
%             elseif vuongRes(i) == 0
%                 annot(i) = "";
%             end
%         elseif strcmp(hypTest,"ad")
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
%             elseif vuongRes(i) == 0
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
% barh(obs,'FaceColor','#d3d3d3');
% hold on
% xline(threshold);
% hold off
% set(gca,'YDir','reverse');
% if fluoOveride
%     ylim([obsLimA+1 obsLimB/2 + 1]);
% else
%     ylim([obsLimA obsLimB]);
% end
% % ylabel('Pressure [dbar]','FontSize',15);
% % set(gca,"YTick",1:1:length(ytix),"YTickLabel",ytix);
% if botCtd == "bot"
%     yticks(0.5:1:19.5);
%     yticklabels(0:10:200);
% elseif botCtd == "ctd"
%     yticks(1:5:55);
%     yticklabels(0:10:150);
% end
% yticklabels({});
% xlabel('No. of Observations',Interpreter='latex',FontSize=13);
% 
% subplot(1,3,[1 2])
% xline(alphaHy,DisplayName='\alpha',LineWidth=1.5,Color="#808080");
% hold on
% if strcmp(hypTest,"ks")
%     if testSel == 4
%         plot(ks(1,:),p,'o-','Color','#a6cee3','DisplayName','Normal','LineWidth',1.5,'MarkerSize',5);
%         plot(ks(2,:),p,'+--','Color','#1f78b4','DisplayName','Lognormal','LineWidth',1.5,'MarkerSize',5);
%         plot(ks(3,:),p,'x-','Color','#b2df8a','DisplayName','Weibull','LineWidth',1.5,'MarkerSize',5);
%         plot(ks(4,:),p,'.--','Color','#33a02c','DisplayName','Gamma','LineWidth',1.5,'MarkerSize',5);
%     elseif testSel == 2
%         plot(ks(1,:),p,'o-','Color','#c51b7d','DisplayName','Normal','LineWidth',1.5,'MarkerSize',5);
%         plot(ks(2,:),p,'+--','Color','#4d9221','DisplayName','Lognormal','LineWidth',1.5,'MarkerSize',5);
%     end
%     xlabel('K-S $p$-value',Interpreter='latex',FontSize=13);
% else
%     if testSel == 4
%         plot(ad(1,:),p,'o-','Color','#a6cee3','DisplayName','Normal','LineWidth',1.5,'MarkerSize',5);
%         plot(ad(2,:),p,'+--','Color','#1f78b4','DisplayName','Lognormal','LineWidth',1.5,'MarkerSize',5);
%         plot(ad(3,:),p,'x-','Color','#b2df8a','DisplayName','Weibull','LineWidth',1.5,'MarkerSize',5);
%         plot(ad(4,:),p,'.--','Color','#33a02c','DisplayName','Gamma','LineWidth',1.5,'MarkerSize',5);
%     elseif testSel == 2
%         for i = 1:n
%             if vuongRes(i) == 1 && ad(1,i) > alphaHy & pV(1,i) > alphaLlr
%                 plot(ad(1,i),p(i),'square','Color','#c51b7d','MarkerSize',15,HandleVisibility='off');
%             elseif vuongRes(i) == 1 && ad(1,i) > alphaHy & pV(1,i) < alphaLlr
%                 plot(ad(1,i),p(i),'square','Color','#c51b7d','MarkerSize',15,'LineWidth',4,HandleVisibility='off');
%             elseif vuongRes(i) == 2 && ad(2,i) > alphaHy & pV(1,i) > alphaLlr
%                 plot(ad(2,i),p(i),'square','Color','#4d9221','MarkerSize',15,HandleVisibility='off');
%             elseif vuongRes(i) == 2 && ad(2,i) > alphaHy & pV(1,i) < alphaLlr
%                 plot(ad(2,i),p(i),'square','Color','#4d9221','MarkerSize',15,'LineWidth',4,HandleVisibility='off');
%             end
%         end
%         plot(ad(1,:),p,'o-','Color','#c51b7d','DisplayName','Normal','LineWidth',1.5,'MarkerSize',5);
%         plot(ad(2,:),p,'+--','Color','#4d9221','DisplayName','Lognormal','LineWidth',1.5,'MarkerSize',5);
%     end
%     xlabel('A-D $p$-value + Vuong Ratio',Interpreter='latex',FontSize=13);
% end
% hold off
% if logAxis == true
%     set(gca, 'XScale', 'log');
%     xline(0.05,':',HandleVisibility='off',LineWidth=1);
%     xline(0.1,':',HandleVisibility='off',LineWidth=1);
% end
% grid minor;
% ylim(limits); xlim([0.6*alphaHy 1]); ylabel('Pressure [dbar]',Interpreter='latex',FontSize=13);
% set(gca,'YDir','reverse');
% legend(Location="best",FontSize=11);
% % yticklabels({});
% % title('K-S p-values');

% subplot(1,5,[3 4])
% zzs = 0.05*ones(n,1);
% for i = 1:n
%     % DEFAULT
%     text(zzs(i),p(i),annot(i),FontSize=11,Color=anClr(i),FontWeight=tmpEmph(i));
% %     % Normal-Lognormal Comparison ONLY
% %     if pV(1,i) > 0.05
% %         text(zzs(i),p(i),annot(i),FontSize=8,Color=anClr(i),FontWeight="bold");
% %     else 
% %         text(zzs(i),p(i),annot(i),FontSize=8,Color=anClr(i));
% %     end
% end
% % % For Normal-Lognormal Comparison ONLY
% % hold on
% % pV(1,obs<threshold) = nan;
% % for i = 1:n
% %     if pV(1,i) > 0.05
% %         scatter(pV(1,i),p(i),[],"black");
% %     end
% % end
% % hold off
% grid minor;
% ylim(limits); set(gca,'YDir','reverse');
% yticklabels({});
% xticklabels({' ' ,' '});
% xlabel('Vuong LLR','FontSize',15);

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

scatter(nan,nan,72,[0.8 0.8 0.8],DisplayName='Data');
hold on
scatter(0,3,72,[0.2 0.2 0.2],'DisplayName','Normal',Marker='pentagram',LineWidth=2.5);
if testSel == 2
    plot(skLogn,kuLogn,'DisplayName','Lognormal','Color','#808080',LineStyle='-',LineWidth=1.3);
    plot(skLognN,kuLognN,'Color','#808080',LineStyle='-',LineWidth=1.3,HandleVisibility='off');
elseif testSel == 4
    plot(skLogn,kuLogn,'DisplayName','Logn.','Color','#1f78b4',LineStyle='--',LineWidth=1.7);
    plot(skLognN,kuLognN,'Color','#1f78b4',LineStyle='--',LineWidth=1.7,HandleVisibility='off');
    plot(skWbl,kuWbl,'DisplayName','Weib.','Color','#b2df8a',LineStyle='-',LineWidth=1.7);
    plot(skWblN,kuWblN,'Color','#b2df8a',LineStyle='-',LineWidth=1.7,HandleVisibility='off');
    plot(skGam,kuGam,'DisplayName','Gam.','Color','#33a02c',LineStyle='--',LineWidth=1.7);
    plot(skGamN,kuGamN,'Color','#33a02c',LineStyle='--',LineWidth=1.7,HandleVisibility='off');
    scatter(2,9,'DisplayName','Exp.',Marker='+',LineWidth=1);
    scatter(0,9/5,'DisplayName','Uni.',Marker='*',LineWidth=1);
    scatter(0,21/5,'DisplayName','Logi.',Marker='.',LineWidth=1);
    scatter(1.1395,5.4,'DisplayName','LEV',Marker='x',LineWidth=1);
end
scatter(sk,ku,72,[0.8 0.8 0.8],HandleVisibility="off");
clr = 1:1:length(p);
scatter(sk,ku,54,clr,"filled","o",HandleVisibility="off");
colormap(gca,cbrewer2("Greens"));
cbar = colorbar;
cbar.Direction = "reverse";
cbar.Ticks = 1:1:length(p);
% cbar.TickLabels = p(1):10:p(end);
cbar.TickLabels = p;
cbar.Label.String = "P [dbar]";
cbar.Label.Position = [0.7 1-0.35];
cbar.Label.Rotation = 0;
% hold on
% add polynomial
[skS,id] = sort(sk);
kuS = ku(id);
[p,S] = polyfit(skS,kuS,2);
[f,delta] = polyval(p,skS,S);
plot(skS,f,'r-',DisplayName="Fit");
plot(skS,f+2*delta,'m--',DisplayName='95% Prediction Interval');
plot(skS,f-2*delta,'m--',HandleVisibility='off');
hold off
grid minor;
ylim([1 kurtLimB]); xlim([skewLimA skewLimB]);
xlabel('Skewness','FontSize',13,'Interpreter','latex'); ylabel('Kurtosis',FontSize=13,Interpreter='latex');
lgd = legend('Location','best');
% title(lgd,'Distributions');
title('L1','Interpreter','latex','FontSize',13);
% sgtitle("L2 chl-$a$ skewness-kurtosis 1988-2021","Interpreter","latex");
exportgraphics(ax,"figures/L1/bottle/log/chla_ad_skKu.png");