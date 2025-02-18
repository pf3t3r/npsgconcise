function [ax,rV,pV,ad,vuongRes] = L2_ctdHelper(X,pIn,maxMld,dcm,testSel,hypTest,limits,threshold,logAxis,legendOn)
% [ax,pL,ks,obs,sk,ku,pV,rV,tr,ad,tr2]
%%L2_helper: this function makes the calculation of KS p-values, skewness,
%%and kurtosis a little more efficient for L2 (sub-mixed layer region that
% is centred on the DCM). 
% INPUTS
% tmp: the .txt datafile with five columns of data (ID, time, time in
% hours, pressure, and variable concentration). Only the ID, pressure, and
% concentration are used.
% maxMld: the maximum mixed layer depth per cruise. Only values shallower
% than this will be considered.
% dcm: pressure of deep chlorophyll maximum (by cruise)
% OUTPUTS
% ax = figure identifier, needed for labelling and export,
% p = pressures where sufficient measurements exist in mixed layer,
% ks = KS test p-values at those pressures,
% obs = observations per depth,
% Sk = skewness at depths where ks is taken,
% Ku = kurtosis at the same depths.

% Set up default values.
if nargin < 10
    legendOn = false;
end
if nargin < 9
    logAxis = true; % true => output p-values in log x-axis, otherwise no log plot.
end
if nargin < 8
    threshold = 50;         % Default threshold = 50 [Mishra et al (2019), 
end                        % Ghasemi & Zahediasl (2012), Ahad et al (2011)]
if nargin < 7
    limits = [-80 140];
end
if nargin < 6
    hypTest = "ks";
end
if nargin < 5
    testSel = 4;
end


n = length(pIn);        % Depth range
nT = length(X(1,:));    % Time range
alphaHy = 0.005;        % Alpha for K-S/A-D p-value
alphaLlr = 0.1;         % Alpha for Vuong LLR p-value


% 1. Extract data beneath ML

pSubml = nan(n,nT);
xSubml = nan(n,nT);
for i = 1:n-1
    %disp(i);
    for j = 1:nT
        if pIn(i) >= maxMld(j)
            pSubml(i,j) = pIn(i+1);
            xSubml(i,j) = X(i+1,j);
        end
    end
end


% 2. Convert p to Lagrangian Coordinates

pL = nan(n,nT);
for i = 1:nT
    pL(:,i) = pSubml(:,i) - dcm(i);
end


% 3. Find KS p-value, Vuong p-value and LLR, skewness, and kurtosis for 
% each depth.

l1 = min(min(pL));
l2 = max(max(pL));

range = l1:2:l2;
rangeLen = 1:1:length(range);
n2 = length(range);
% disp(n2);

% Automate limit setting
bottom = floor((l1 - limits(1))./2);
top = floor((l2 - limits(2))./2); 
% indexOne = find(pL==limits(1));
% indexTwo = find(pL==limits(2));
% disp(bottom);
% disp(top);

ks = nan(5,n2);
ad = nan(4,n2);
rV = nan(10,n2); 
pV = nan(10,n2);
obs = nan(1,n2);
sk = nan(1,n2);
ku = nan(1,n2);
%adH = nan(1,n2); adP = nan(1,n2); % optional parameters for A-D Test

for i = rangeLen
    tmp = xSubml(pL==range(i));
    %disp(i);
    tmp(tmp<=0) = nan;
    tmp(isnan(tmp)) = [];
    obs(i) = length(tmp);
    if length(tmp) > 3
        gammaParams = mle(tmp,"distribution","Gamma");
        pdG = makedist("Gamma",gammaParams(1),gammaParams(2));
        [~,ks(:,i),~] = statsplot2(tmp,'noplot');
        [~,ad(2,i)] = adtest(tmp,"Distribution","logn");
        [~,ad(1,i)] = adtest(tmp,"Distribution","norm");
        [~,ad(3,i)] = adtest(tmp,"Distribution","weibull");
        [~,ad(4,i)] = adtest(tmp,Distribution=pdG,MCTol=0.05);
        [rV(:,i),pV(:,i)] = bbvuong(tmp);
        %[adH(i),adP(i)] = lillietest(log(tmp),MCTol=1e-2); % optional Lilliefors Test
        sk(i) = skewness(tmp);
        ku(i) = kurtosis(tmp);
    end
end

% Set results = nan if they do not meet the threshold
for i = rangeLen
    if obs(i) < threshold
        ks(:,i) = nan;
        ad(:,i) = nan;
        ku(i) = nan;
        pV(:,i) = nan;
        rV(:,i) = nan;
        sk(i) = nan;
    end
end


% 4. Compare all Vuong Test LLR results to give a single best distribution
% for each depth. Comment / uncomment only one of 4.a or 4.b.

vuongRes = zeros(1,n2);
annot = strings(1,n2);
anClr = strings(1,n2);
anClr(cellfun(@isempty,anClr)) = '#FFFFFF';
tmpEmph = strings(1,n2); tmpEmph(cellfun(@isempty,tmpEmph)) = 'bold';
rV(isnan(rV)) = 0;

% 4.a. Vuong: Normal vs Lognormal vs Weibull vs Gamma
if testSel==4 && hypTest == "ks"
    for i = 1:n2
        if rV(1,i) > 0  & rV(2,i) > 0 & rV(3,i) > 0 & ks(1,i) > alphaHy
            if length(find(ks(:,i)>alphaHy)) == 1
                tmp = "";
            else
                vuongRes(i) = 1;
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
        elseif rV(1,i) < 0 & rV(5,i) > 0 & rV(6,i) > 0 & ks(2,i) > alphaHy
            if length(find(ks(:,i)>alphaHy)) == 1
                tmp = "";
            else
                vuongRes(i) = 2;
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
        elseif rV(2,i) < 0 & rV(5,i) < 0 & rV(8,i) > 0 & ks(3,i) > alphaHy
            if length(find(ks(:,i)>alphaHy)) == 1
                tmp = "";
            else
                vuongRes(i) = 3;
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
        elseif rV(3,i) < 0 & rV(6,i) < 0 & rV(8,i) < 0 & ks(4,i) > alphaHy
            if length(find(ks(:,i)>alphaHy)) == 1
                tmp = "";
            else
                vuongRes(i) = 4;
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
        end
    end
    rV(rV==0) = nan;
elseif testSel == 4 && hypTest == "ad"
    %%%
    for i = 1:n2
        if rV(1,i) > 0 & rV(2,i) > 0 & rV(3,i) > 0 & ad(1,i) > alphaHy
            if length(find(ad(:,i)>alphaHy)) == 1
                tmp = "";
            else
                vuongRes(i) = 1;
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
        elseif rV(1,i) < 0 & rV(5,i) > 0 & rV(6,i) > 0 & ad(2,i) > alphaHy
            if length(find(ad(:,i)>alphaHy)) == 1
                tmp = "";
            else
                vuongRes(i) = 2;
                tmp = "Lognormal";
            end
            anClr(i) = '#1f78b4';
            if pV(1,i) > alphaLlr & ad(1,i) > alphaHy
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
        elseif rV(2,i) < 0 & rV(5,i) < 0 & rV(8,i) > 0 & ad(3,i) > alphaHy
            if length(find(ad(:,i)>alphaHy)) == 1
                tmp = "";
            else
                vuongRes(i) = 3;
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
        elseif rV(3,i) < 0 & rV(6,i) < 0 & rV(8,i) < 0 & ad(4,i) > alphaHy
            if length(find(ad(:,i)>alphaHy)) == 1
                tmp = "";
            else
                vuongRes(i) = 4;
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
        end
    end
    rV(rV==0) = nan;
    %%%
elseif testSel == 2
    % 4.b. Vuong: Normal Vs. Lognormal Only
    for i = 1:n2
        if rV(1,i) > 0 
            vuongRes(i) = 1;
            annot(i) = "Normal";
            anClr(i) = '#a6cee3';
            if pV(1,i) > alphaLlr
                tmpEmph(i) = 'normal';
            end
        elseif rV(1,i) < 0
            vuongRes(i) = 2;
            annot(i) = "Lognormal";
            anClr(i) = '#1f78b4';
            if pV(1,i) > alphaLlr
                tmpEmph(i) = 'normal';
            end
        else
            annot(i) = "";
        end
    end
    rV(rV==0) = nan;
end

% % 5. Generate theoretical skewness-kurtosis curves for the...
% % Lognormal family,
% sigTh = linspace(0,1,1000);
% for i = 1:length(sigTh)
%     skLogn(i) = (exp(sigTh(i)^2) + 2)*(sqrt(exp(sigTh(i)^2) - 1));
%     kuLogn(i) = exp(4*sigTh(i)^2) + 2*exp(3*sigTh(i)^2) + 3*exp(2*sigTh(i)^2) - 3;
% end
% 
% % Gamma family, and
% kTh = linspace(0.04,3000,1500000);
% for i = 1:length(kTh)
%     skGam(i) = 2/sqrt(kTh(i));
%     kuGam(i) = 6/kTh(i) + 3;
% end
% 
% % Weibull family
% kWbl = linspace(0.1,3.5,10000);
% for i = 1:length(kWbl)
%     skWbl(i) = ( gamma(1 + 3/kWbl(i)) - 3*gamma(1 + 1/kWbl(i))*gamma(1 + 2/kWbl(i)) + 2*(gamma(1 + 1/kWbl(i)))^3 ) ./ ...
%         ( gamma(1 + 2/kWbl(i)) -  (gamma(1 + 1/kWbl(i)))^2 )^(3/2);
%     kuWbl(i) = ( gamma(1 + 4/kWbl(i)) - 4*gamma(1 + 1/kWbl(i))*gamma(1 + 3/kWbl(i)) + 6*( (gamma(1 + 1/kWbl(i)) )^2)*gamma(1 + 2/kWbl(i)) - 3*( (gamma(1 + 1/kWbl(i)))^4 ) ) ./ ...
%        ( gamma(1 + 2/kWbl(i)) - ( gamma(1 + 1/kWbl(i)) )^2 )^2;
% end

% % Negative Distributions
% skLognN = -skLogn;
% kuLognN = kuLogn;
% skGamN= -skGam;
% kuGamN = kuGam;
% skWblN = -skWbl;
% kuWblN = kuWbl;

% 6. Plot everything.

ax = figure;

subplot(1,3,3)
barh(obs(rangeLen(1):rangeLen(end)),'FaceColor','#d3d3d3');
hold on
xline(threshold);
hold off
set(gca,'YDir','reverse');
% set(gca,'XDir','reverse');
% ylabel('Pressure [dbar]',FontSize=13);
xlabel('No. of Obs.','FontSize',13,Interpreter='latex');
set(gca,"YTick",2:5:n2,"YTickLabel",range(2):10:range(end));
ylim([rangeLen(1-bottom) rangeLen(end-top)]);
yticklabels({})
% title('No. of Observations');

subplot(1,3,[1 2])
xline(alphaHy,'-','\alpha',LineWidth=1.5,Color="#808080",HandleVisibility="off");
hold on
if strcmp(hypTest,"ks")
    if testSel==4
        plot(ks(1,:),range,'o-','Color','#a6cee3','DisplayName','Normal','LineWidth',1.5,'MarkerSize',5);
        plot(ks(2,:),range,'+--','Color','#1f78b4','DisplayName','Lognormal','LineWidth',1.5,'MarkerSize',5);
        plot(ks(3,:),range,'x-','Color','#b2df8a','DisplayName','Weibull','LineWidth',1.5,'MarkerSize',5);
        plot(ks(4,:),range,'.--','Color','#33a02c','DisplayName','Gamma','LineWidth',1.5,'MarkerSize',5);
    elseif testSel==2
        plot(ks(1,:),range,'o-','Color','#a6cee3','DisplayName','Normal','LineWidth',1.5,'MarkerSize',5);
        plot(ks(2,:),range,'+--','Color','#1f78b4','DisplayName','Lognormal','LineWidth',1.5,'MarkerSize',5);
    end
    xlabel('K-S $p$-value','Interpreter','latex',FontSize=15);
elseif strcmp(hypTest,"ad")
    if testSel==4
        plot(ad(1,:),range,'o-','Color','#a6cee3','DisplayName','Normal','LineWidth',1.5,'MarkerSize',5);
        plot(ad(2,:),range,'+--','Color','#1f78b4','DisplayName','Lognormal','LineWidth',1.5,'MarkerSize',5);
        plot(ad(3,:),range,'x-','Color','#b2df8a','DisplayName','Weibull','LineWidth',1.5,'MarkerSize',5);
        plot(ad(4,:),range,'.--','Color','#33a02c','DisplayName','Gamma','LineWidth',1.5,'MarkerSize',5);
    elseif testSel==2
        for i = 1:n
            if vuongRes(i) == 1 && ad(1,i) > alphaHy & pV(1,i) > alphaLlr && ad(2,i) > alphaHy
                plot(ad(1,i),range(i),'square','Color','#c51b7d','MarkerSize',15,HandleVisibility='off');
            elseif vuongRes(i) == 1 && ad(1,i) > alphaHy & pV(1,i) < alphaLlr && ad(2,i) > alphaHy
                plot(ad(1,i),range(i),'square','Color','#c51b7d','MarkerSize',15,'LineWidth',4,HandleVisibility='off');
            elseif vuongRes(i) == 2 && ad(2,i) > alphaHy & pV(1,i) > alphaLlr && ad(1,i) > alphaHy
                plot(ad(2,i),range(i),'square','Color','#4d9221','MarkerSize',15,HandleVisibility='off');
            elseif vuongRes(i) == 2 && ad(2,i) > alphaHy & pV(1,i) < alphaLlr && ad(1,i) > alphaHy
                plot(ad(2,i),range(i),'square','Color','#4d9221','MarkerSize',15,'LineWidth',4,HandleVisibility='off');
            end
        end
        plot(ad(1,:),range,'o-','Color','#c51b7d','DisplayName','Normal','LineWidth',1.5,'MarkerSize',5);
        plot(ad(2,:),range,'o-','Color','#4d9221','DisplayName','Lognormal','LineWidth',1.5,'MarkerSize',5);
    end
    xlabel('A-D $p$-value',FontSize=13,Interpreter='latex');
end
hold off
grid minor;
if logAxis == true
    set(gca, 'XScale', 'log');
    %xline(0.05,'--',HandleVisibility='off');
    %xline(0.1,'--',HandleVisibility='off');
end
ylim(limits); xlim([0.1*alphaHy 1]); 
set(gca,"YTick",limits(1):10:limits(2),"YTickLabel",limits(1):10:limits(2));
% yticklabels({});
ylabel("Pressure [dbar]",Interpreter="latex",FontSize=13);
set(gca,'YDir','reverse');
if legendOn == true
    legend(Location="west");
end
% title('K-S Test');

% zzs = 0.25*ones(n2,1);
% subplot(1,6,4)
% for i = 1-bottom:n2-top
%     text(zzs(i),range(i),annot(i),FontSize=8,Color=anClr(i),FontWeight=tmpEmph(i));
% end
% % ylim([l1+40 l2-60]); 
% ylim(limits);
% set(gca,"YTick",limits(1):10:limits(2),"YTickLabel",limits(1):10:limits(2));
% % yticklabels({});
% set(gca,'YDir','reverse');
% xticklabels({' ' ,' '});
% xlabel('Vuong LLR','FontSize',15);
% % set(gca,'xtick',[]);
% % title('Vuong LLR');
% 
% tmp = [];
% for i = 1:n2
%     if ~isnan(sum(ks(:,i)))
%         tmp = [tmp i];
%     end
% end
% % disp(tmp);
% tr2 = range(tmp);
% sk2 = sk(tmp);
% ku2 = ku(tmp);
% clear tmp;
% 
% if ismember(tr2,limits(1))
%     kLimA = find(tr2==limits(2));
% else
%     kLimA = 1;
% end
% 
% kLim = [kLimA find(tr2==limits(2))];
% 
% tr2 = tr2(kLim(1):kLim(2));
% sk2 = sk2(kLim(1):kLim(2));
% ku2 = ku2(kLim(1):kLim(2));
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
% % error bars for sk-ku
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
% subplot(1,6,[5 6])
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
% scatter(sk2,ku2,[],[0.8 0.8 0.8],HandleVisibility="off");
% clr = 1:1:length(tr2);
% scatter(sk2,ku2,24,clr,"filled","o",HandleVisibility="off");
% colormap(gca,cbrewer2("RdYlBu"));
% cbar = colorbar;
% cbar.Direction = "reverse";
% cbar.Ticks = 1:5:(kLim(2)-kLim(1)+5);
% cbar.TickLabels = limits(1):10:limits(2);
% % cbar.Ticks = 1:5:(limits(2)-limits(1))/2;
% % cbar.TickLabels = limits(1):10:limits(2);
% % cbar.Ticks = 3:10:length(tr2);
% % cbar.TickLabels = tr2(3):20:tr2(end);
% cbar.Label.String = "P [dbar]";
% hold off
% grid minor;
% ylim([1 kurtLimB]); xlim([skewLimA skewLimB]);
% xlabel('Skewness','FontSize',15); ylabel('Kurtosis',FontSize=15);
% lgd = legend('Location','best');
% title(lgd,'Distributions');
% % title('Skewness vs. Kurtosis');
% 
% sk = sk2; ku = ku2; tr = tr2;

end