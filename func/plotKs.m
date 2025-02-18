function [] = plotKs(tr,ks,obs,sk,ku,obsLimA,obsLimB,EulLan,threshold,vuongRes,pV,limitOveride,fluoOveride,hypTest,ad,testSel,botCtd,logAxis)
%plotKs
% INPUT: 
% OUTPUT: 

if nargin < 18
    logAxis = true;
end

if nargin < 17
    botCtd = "bot";
end

if nargin < 16
    testSel = 4;
end

if nargin < 14
    hypTest = "ks";
end

if nargin < 13
    fluoOveride = false;
end

% 'unc' unused: add to func input if used in future
% if nargin < 13
%     unc = nan(70,16);
% end

if nargin < 10
    vuongRes = 0;
end

if nargin < 9
    threshold = 50;
end

if nargin < 8
    EulLan = true;
end

if nargin < 6
    obsLimA = 0.5;
    %obsLimB = length(tr);
    obsLimB = length(obs) + 0.5;
end

if EulLan
    limits = [0 200];
    if obsLimB-obsLimA > 20
        ytix = 5:5:205;
    else
        ytix = 5:10:205;
    end
else
    limits = [-150 100];
    ytix = tr;
end

if nargin >= 12
    limits = limitOveride;
end

if fluoOveride
    ytix = obsLimA:2:obsLimB;
end

n = length(tr);

% Create Annotations for Vuong's Test Results
annot = strings(1,n);
anClr = strings(1,n);
anClr(cellfun(@isempty,anClr)) = '#FFFFFF';
tmpEmph = strings(1,n); tmpEmph(cellfun(@isempty,tmpEmph)) = 'bold';
% Default Case

alphaHy = 0.005;
alphaLlr = 0.10;

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
            % A-D
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
    % Normal-Lognormal Case
    for i = 1:n
        if strcmp(hypTest,"ks")
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
            elseif vuongRes(i) == 0
                annot(i) = "";
            end
        elseif strcmp(hypTest,"ad")
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
            elseif vuongRes(i) == 0
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

% % Beta family: generate theoretical skewness and kurtosis
% muB = linspace(0.5,5,100);
% nuB = linspace(0.5,5,100);
% alphaB = muB.*nuB; betaB = (1 - muB).*nuB;
% for i = 1:length(betaB)
%     skBet(i) = 2*(betaB(i) - alphaB(i)) * sqrt(alphaB(i) + betaB(i) + 1) ./ ...
%         sqrt(betaB(i)*alphaB(i)) * (betaB(i) + alphaB(i) + 2);
%     kuBet(i) = 3*(betaB(i) + alphaB(i) + 1) * ( 2*(betaB(i) + alphaB(i))^2 + alphaB(i)*betaB(i)*(betaB(i) + alphaB(i) - 6) ) ./ ...
%         betaB(i)*alphaB(i)*(betaB(i) + alphaB(i) + 2)*(betaB(i) + alphaB(i) + 3);
% end

% Loglogistic family: generate theoretical skewness and kurtosis
% Nothings shows up: Maybe I need to tune the shape parameter a bit (??)
% cLgi = linspace(2,4,10000);
% for i = 1:length(cLgi)
%     skLgi(i) = ( 2*(pi^2)*(csc(pi/cLgi(i)))^3 - 6*cLgi(i)*pi*csc(pi/cLgi(i))*csc(2*pi/cLgi(i)) + 3*(cLgi(i)^2)*csc(3*pi/cLgi(i)) ) ./ ...
%         sqrt(pi)*((-pi*(csc(pi/cLgi(i))^2)) + 2*cLgi(i)*csc(2*pi/cLgi(i)))^(3/2);
%     kuLgi(i) = ( -3*pi^2*(csc(pi/cLgi(i))^4) - 12*pi*cLgi(i)^2*csc(pi/cLgi(i))*csc(3*pi/cLgi(i)) + 4*cLgi(i)^3*csc(4*pi/cLgi(i)) + 6*cLgi(i)*pi^2*(csc(pi/cLgi(i))^3)*sec(pi/cLgi(i)) ) ./ ...
%         pi*( -pi*(csc(pi/cLgi(i))^2) + 2*cLgi(i)*csc(2*pi/cLgi(i)) )^2;
% end

% Negative Distributions
skLognN = -skLogn;
kuLognN = kuLogn;
% skGamN= -skGam;
% kuGamN = kuGam;
% skWblN = -skWbl;
% kuWblN = kuWbl;

subplot(1,3,3)
barh(obs,'FaceColor','#d3d3d3');
hold on
xline(threshold);
hold off
set(gca,'YDir','reverse');
if fluoOveride
    ylim([obsLimA+1 obsLimB/2 + 1]);
else
    ylim([obsLimA obsLimB]);
end
% ylabel('Pressure [dbar]','FontSize',15);
% set(gca,"YTick",1:1:length(ytix),"YTickLabel",ytix);
if botCtd == "bot"
    yticks(0.5:1:19.5);
    yticklabels(0:10:200);
elseif botCtd == "ctd"
    yticks(1:5:55);
    yticklabels(0:10:150);
end
yticklabels({});
xlabel('No. of Obs.',Interpreter='latex',FontSize=13);

subplot(1,3,[1 2])
xline(alphaHy,'-','\color{black}\alpha=0.005',LineWidth=1.5,Color="#808080",HandleVisibility="off",LabelOrientation="horizontal",LabelHorizontalAlignment="center",FontSize=13);
hold on
if strcmp(hypTest,"ks")
    if testSel == 4
        plot(ks(1,:),tr,'o-','Color','#a6cee3','DisplayName','Normal','LineWidth',1.5,'MarkerSize',5);
        plot(ks(2,:),tr,'+-','Color','#1f78b4','DisplayName','Lognormal','LineWidth',1.5,'MarkerSize',5);
        plot(ks(3,:),tr,'x-','Color','#b2df8a','DisplayName','Weibull','LineWidth',1.5,'MarkerSize',5);
        plot(ks(4,:),tr,'.-','Color','#33a02c','DisplayName','Gamma','LineWidth',1.5,'MarkerSize',5);
    elseif testSel == 2
        plot(ks(1,:),tr,'o-','Color','#c51b7d','DisplayName','Normal','LineWidth',1.5,'MarkerSize',5);
        plot(ks(2,:),tr,'+--','Color','#4d9221','DisplayName','Lognormal','LineWidth',1.5,'MarkerSize',5);
    end
    xlabel('K-S $p$-value',Interpreter='latex',FontSize=13);
else
    if testSel == 4
        plot(ad(1,:),tr,'o-','Color','#c51b7d','DisplayName','Normal','LineWidth',1.5,'MarkerSize',5);
        plot(ad(2,:),tr,'+-','Color','#4d9221','DisplayName','Lognormal','LineWidth',1.5,'MarkerSize',5);
        plot(ad(3,:),tr,'x-','Color','#d3d3d3','DisplayName','Weibull','LineWidth',1.5,'MarkerSize',5);
        plot(ad(4,:),tr,'.-','Color','#000000','DisplayName','Gamma','LineWidth',1.5,'MarkerSize',5);
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
        plot(nan,nan,'square','Color','#808080','MarkerSize',15,'DisplayName','V-LLR best fit (p > 0.1)');        
        plot(nan,nan,'square','Color','#808080','MarkerSize',15,'LineWidth',4,'DisplayName','V-LLR best fit (p < 0.1)');        
        plot(ad(1,:),tr,'o-','Color','#c51b7d','DisplayName','Normal','LineWidth',1.5,'MarkerSize',5);
        plot(ad(2,:),tr,'o-','Color','#4d9221','DisplayName','Lognormal','LineWidth',1.5,'MarkerSize',5);
    end
    xlabel('A-D $p$-value',Interpreter='latex',FontSize=13);
end
hold off
if logAxis == true
    set(gca, 'XScale', 'log');
    %xline(0.05,':',HandleVisibility='off',LineWidth=1);
    %xline(0.1,':',HandleVisibility='off',LineWidth=1);
end
grid minor;
ylim(limits); xlim([0.1*alphaHy 1]); ylabel('Pressure [dbar]',Interpreter='latex',FontSize=13);
set(gca,'YDir','reverse');
legend('Position',[0.4 0.7 0.07 0.12],FontSize=11);
% yticklabels({});
% title('K-S p-values');

% subplot(1,5,[3 4])
% zzs = 0.05*ones(n,1);
% for i = 1:n
%     % DEFAULT
%     text(zzs(i),tr(i),annot(i),FontSize=11,Color=anClr(i),FontWeight=tmpEmph(i));
% %     % Normal-Lognormal Comparison ONLY
% %     if pV(1,i) > 0.05
% %         text(zzs(i),tr(i),annot(i),FontSize=8,Color=anClr(i),FontWeight="bold");
% %     else 
% %         text(zzs(i),tr(i),annot(i),FontSize=8,Color=anClr(i));
% %     end
% end
% % % For Normal-Lognormal Comparison ONLY
% % hold on
% % pV(1,obs<threshold) = nan;
% % for i = 1:n
% %     if pV(1,i) > 0.05
% %         scatter(pV(1,i),tr(i),[],"black");
% %     end
% % end
% % hold off
% grid minor;
% ylim(limits); set(gca,'YDir','reverse');
% yticklabels({});
% xticklabels({' ' ,' '});
% xlabel('Vuong LLR','FontSize',15);

% kurtLimB = 10; skewLimA = 0; skewLimB = 2.5;
% if max(ku) > 10 & min(sk) < 0
%     kurtLimB = max(ku) + 1;
%     skewLimA = min(sk) - 0.1;
%     skewLimB = max(sk) + 0.1;
% elseif max(ku) > 10
%     kurtLimB = max(ku) + 1;
%     skewLimB = max(sk) + 0.1;
% elseif min(sk) < 0 
%     skewLimA = min(sk) - 0.1;
% elseif max(sk) > 2.5
%     kurtLimB = max(ku) + 1;
%     skewLimB = max(sk) + 0.1;
% else 
%     kurtLimB = 10;
%     skewLimA = 0;
%     skewLimB = 2.5;
% end
% 
% % % error bars for sk-ku
% % obsTmp = obs(obs>=50);
% % disp(length(obsTmp));
% % disp(n);
% % ynegL = nan(n,1); yposL = nan(n,1); xnegL = nan(n,1); xposL = nan(n,1);
% % ynegN = nan(n,1); yposN = nan(n,1); xnegN = nan(n,1); xposN = nan(n,1);
% % ynegG = nan(n,1); yposG = nan(n,1); xnegG = nan(n,1); xposG = nan(n,1);
% % ynegW = nan(n,1); yposW = nan(n,1); xnegW = nan(n,1); xposW = nan(n,1);
% % for i = 1:n
% %     if obsTmp(i) > 300
% %         ynegL(i) = unc(60,5); yposL(i) = unc(60,6);
% %         %ynegL(i) = -1.02; yposL(i) = 0.92;
% %         xnegL(i) = unc(60,7); xposL(i) = unc(60,8);
% %         %xnegL(i) = -0.23; xposL(i) = 0.22;
% %         ynegN(i) = unc(60,1); yposN(i) = unc(60,2);
% %         xnegN(i) = unc(60,3); xposN(i) = unc(60,4);
% %         ynegG(i) = unc(60,9); yposG(i) = unc(60,10);
% %         xnegG(i) = unc(60,11); xposG(i) = unc(60,12);
% %         ynegW(i) = unc(60,13); yposW(i) = unc(60,14);
% %         xnegW(i) = unc(60,15); xposW(i) = unc(60,16);
% %     elseif obsTmp(i) > 250
% %         ynegL(i) = unc(50,5); yposL(i) = unc(50,6);
% %         %ynegL(i) = -1.06; yposL(i) = 0.96;
% %         xnegL(i) = unc(50,7); xposL(i) = unc(50,8);
% %         %xnegL(i) = -0.24; xposL(i) = 0.25;
% %         ynegN(i) = unc(50,1); yposN(i) = unc(50,2);
% %         xnegN(i) = unc(50,3); xposN(i) = unc(50,4);
% %         ynegG(i) = unc(50,9); yposG(i) = unc(50,10);
% %         xnegG(i) = unc(50,11); xposG(i) = unc(50,12);
% %         ynegW(i) = unc(50,13); yposW(i) = unc(50,14);
% %         xnegW(i) = unc(50,15); xposW(i) = unc(50,16);
% %     elseif obsTmp(i) > 200
% %         ynegL(i) = unc(40,5); yposL(i) = unc(40,6);
% %         %ynegL(i) = -1.13; yposL(i) = 0.96;
% %         xnegL(i) = unc(40,7); xposL(i) = unc(40,8);
% %         %xnegL(i) = -0.27; xposL(i) = 0.26;
% %         ynegN(i) = unc(40,1); yposN(i) = unc(40,2);
% %         xnegN(i) = unc(40,3); xposN(i) = unc(40,4);
% %         ynegG(i) = unc(40,9); yposG(i) = unc(40,10);
% %         xnegG(i) = unc(40,11); xposG(i) = unc(40,12);
% %         ynegW(i) = unc(40,13); yposW(i) = unc(40,14);
% %         xnegW(i) = unc(40,15); xposW(i) = unc(40,16);
% %     elseif obsTmp(i) > 150
% %         ynegL(i) = unc(30,5); yposL(i) = unc(30,6);
% %         %ynegL(i) = -1.22; yposL(i) = 1.08;
% %         xnegL(i) = unc(30,7); xposL(i) = unc(30,8);
% %         %xnegL(i) = -0.30; xposL(i) = 0.29;
% %         ynegN(i) = unc(30,1); yposN(i) = unc(30,2);
% %         xnegN(i) = unc(30,3); xposN(i) = unc(30,4);
% %         ynegG(i) = unc(30,9); yposG(i) = unc(30,10);
% %         xnegG(i) = unc(30,11); xposG(i) = unc(30,12);
% %         ynegW(i) = unc(30,13); yposW(i) = unc(30,14);
% %         xnegW(i) = unc(30,15); xposW(i) = unc(30,16);
% %     elseif obsTmp(i) > 100
% %         ynegL(i) = unc(20,5); yposL(i) = unc(20,6);
% %         %ynegL(i) = -1.30; yposL(i) = 1.15;
% %         xnegL(i) = unc(20,7); xposL(i) = unc(20,8);
% %         %xnegL(i) = -0.34; xposL(i) = 0.33;
% %         ynegN(i) = unc(20,1); yposN(i) = unc(20,2);
% %         xnegN(i) = unc(20,3); xposN(i) = unc(20,4);
% %         ynegG(i) = unc(20,9); yposG(i) = unc(20,10);
% %         xnegG(i) = unc(20,11); xposG(i) = unc(20,12);
% %         ynegW(i) = unc(20,13); yposW(i) = unc(20,14);
% %         xnegW(i) = unc(20,15); xposW(i) = unc(20,16);
% %     else
% %         ynegL(i) = unc(10,5); yposL(i) = unc(10,6);
% %         %ynegL(i) = -1.24; yposL(i) = 1.11;
% %         xnegL(i) = unc(10,7); xposL(i) = unc(10,8);
% %         %xnegL(i) = -0.51; xposL(i) = 0.38;
% %         ynegN(i) = unc(10,1); yposN(i) = unc(10,2);
% %         xnegN(i) = unc(10,3); xposN(i) = unc(10,4);
% %         ynegG(i) = unc(10,9); yposG(i) = unc(10,10);
% %         xnegG(i) = unc(10,11); xposG(i) = unc(10,12);
% %         ynegW(i) = unc(10,13); yposW(i) = unc(10,14);
% %         xnegW(i) = unc(10,15); xposW(i) = unc(10,16);
% %     end
% % end
% 
% subplot(1,6,[4.1 5])
% scatter(0,3,[],[0.7725	0.1059	0.4902],'DisplayName','Norm.',Marker='o',LineWidth=3);
% hold on
% if testSel == 2
%     plot(skLogn,kuLogn,'DisplayName','Logn.','Color','#4d9221',LineStyle='--',LineWidth=1.7);
%     plot(skLognN,kuLognN,'Color','#4d9221',LineStyle='--',LineWidth=1.7,HandleVisibility='off');
% elseif testSel == 4
%     plot(skLogn,kuLogn,'DisplayName','Logn.','Color','#1f78b4',LineStyle='--',LineWidth=1.7);
%     plot(skLognN,kuLognN,'Color','#1f78b4',LineStyle='--',LineWidth=1.7,HandleVisibility='off');
%     plot(skWbl,kuWbl,'DisplayName','Weib.','Color','#b2df8a',LineStyle='-',LineWidth=1.7);
%     plot(skWblN,kuWblN,'Color','#b2df8a',LineStyle='-',LineWidth=1.7,HandleVisibility='off');
%     plot(skGam,kuGam,'DisplayName','Gam.','Color','#33a02c',LineStyle='--',LineWidth=1.7);
%     plot(skGamN,kuGamN,'Color','#33a02c',LineStyle='--',LineWidth=1.7,HandleVisibility='off');
%     scatter(2,9,'DisplayName','Exp.',Marker='+',LineWidth=1);
%     scatter(0,9/5,'DisplayName','Uni.',Marker='*',LineWidth=1);
%     scatter(0,21/5,'DisplayName','Logi.',Marker='.',LineWidth=1);
%     scatter(1.1395,5.4,'DisplayName','LEV',Marker='x',LineWidth=1);
% end
% % errorbar(sk,ku,ynegL,yposL,xnegL,xposL,'o','Color','#1f78b4','HandleVisibility','off',LineWidth=1.9);
% % errorbar(sk,ku,ynegG,yposG,xnegG,xposG,'o','Color','#33a02c','HandleVisibility','off',LineWidth=1.6);
% % errorbar(sk,ku,ynegW,yposW,xnegW,xposW,'o','Color','#b2df8a','HandleVisibility','off',LineWidth=1.3);
% % errorbar(sk,ku,ynegN,yposN,xnegN,xposN,'o','Color',[0.6509803921568628 0.807843137254902 0.8901960784313725],'HandleVisibility','off',LineWidth=1);
% scatter(sk,ku,[],[0.8 0.8 0.8],HandleVisibility="off");
% clr = 1:1:length(tr);
% scatter(sk,ku,24,clr,"filled","o",HandleVisibility="off");
% colormap(gca,cbrewer2("PiYG"));
% cbar = colorbar;
% cbar.Direction = "reverse";
% cbar.Ticks = 1:1:length(tr);
% % cbar.TickLabels = tr(1):10:tr(end);
% cbar.TickLabels = tr;
% cbar.Label.String = "P [dbar]";
% % hold on
% hold off
% grid minor;
% ylim([1 kurtLimB]); xlim([skewLimA skewLimB]);
% xlabel('Skewness','FontSize',15); ylabel('Kurtosis',FontSize=15);
% lgd = legend('Location','best');
% title(lgd,'Distributions');
% % title('Skewness vs. Kurtosis');

end