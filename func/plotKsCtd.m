function plotKsCtd(seasYrly,ks1,ks2,ks3,ks4,p,MLD,dcmMu,dcmSigma,dcmIsoMu,dcmIsoSigma)
%plotKsCtd Plot quickly the results of the Kolmogorov-Smirnov test for four
% possible situations: These can either be (i) Eulerian Isopycnal, (ii) 
% Eulerian Pressure, (iii) Lagrangian Isopycnal, (iv) Lagrangian Pressure.
% OR: (i) Winter, (ii) Spring, (iii) Summer, (iv) Autumn. This function 
% will be updated in the future to account for KS tests run in actual 
% isopycnal coordinates, and we will rename the current isopycnal test.

% Smooth KS Results over 10 dbar
ks1 = movmean(ks1,5,2);
ks2 = movmean(ks2,5,2);
ks3 = movmean(ks3,5,2);
ks4 = movmean(ks4,5,2);

if seasYrly == "yrly"
    % Isopycnal Eulerian
    subplot(1,4,1)
    plot(ks1(1,:),p,'o-','Color','#a6cee3','DisplayName','Normal','LineWidth',1.5,'MarkerSize',5);
    hold on
    plot(ks1(2,:),p,'+--','Color','#1f78b4','DisplayName','Lognormal','LineWidth',1.5,'MarkerSize',5);
    plot(ks1(3,:),p,'x-','Color','#b2df8a','DisplayName','Weibull','LineWidth',1.5,'MarkerSize',5);
    plot(ks1(4,:),p,'.--','Color','#33a02c','DisplayName','Gamma','LineWidth',1.5,'MarkerSize',5);
    if nargin==11
        yline(dcmIsoMu,'Color',[0.3 0.3 0.3],'LineWidth',5,'alpha',0.5,'DisplayName','Mean DCM','HandleVisibility','off');
        yline(dcmIsoMu+dcmIsoSigma,'k:','LineWidth',2,'DisplayName','SD','HandleVisibility','off');
        yline(dcmIsoMu-dcmIsoSigma,'k:','LineWidth',2,'HandleVisibility','off');
    end
    hold off
    set(gca,'YDir','reverse');
    h_leg = legend('Location','best');
    ylim([0 250]);
    xlabel('p-value');
    ylabel('Pressure [db]');
    title('Isopycnal Eulerian');
    
    % Eulerian
    subplot(1,4,2)
    plot(ks2(1,:),p,'o-','Color','#a6cee3','DisplayName','Normal','LineWidth',1.5,'MarkerSize',5);
    hold on
    plot(ks2(2,:),p,'+--','Color','#1f78b4','DisplayName','Lognormal','LineWidth',1.5,'MarkerSize',5);
    plot(ks2(3,:),p,'x-','Color','#b2df8a','DisplayName','Weibull','LineWidth',1.5,'MarkerSize',5);
    plot(ks2(4,:),p,'.--','Color','#33a02c','DisplayName','Gamma','LineWidth',1.5,'MarkerSize',5);
    yline(mean(MLD),'Color',[0.4 0.4 0.4],'LineWidth',3,'DisplayName','Average MLD','HandleVisibility','off');
    if nargin>7
        yline(dcmMu,'Color',[0.3 0.3 0.3],'LineWidth',5,'DisplayName','Mean DCM','HandleVisibility','off');
    end
    if nargin>8
        yline(dcmMu+dcmSigma,'k:','LineWidth',2,'DisplayName','SD','HandleVisibility','off');
        yline(dcmMu-dcmSigma,'k:','LineWidth',2,'DisplayName','SD','HandleVisibility','off');
    end
    hold off
    set(gca,'YDir','reverse');
    ylim([0 250]);
    xlabel('p-value');
    title('Eulerian');
    
    % Isopycnal Lagrangian
    subplot(1,4,3)
    plot(ks3(1,:),p-129,'o-','Color','#a6cee3','DisplayName','Normal','LineWidth',1.5,'MarkerSize',5);
    hold on
    plot(ks3(2,:),p-129,'+--','Color','#1f78b4','DisplayName','Lognormal','LineWidth',1.5,'MarkerSize',5);
    plot(ks3(3,:),p-129,'x-','Color','#b2df8a','DisplayName','Weibull','LineWidth',1.5,'MarkerSize',5);
    plot(ks3(4,:),p-129,'.--','Color','#33a02c','DisplayName','Gamma','LineWidth',1.5,'MarkerSize',5);
    hold off
    set(gca,'YDir','reverse');
    xlabel('p-value');
    ylim([-125 125]);
    title('Isopycnal Lagrangian');
    
    % Lagrangian
    subplot(1,4,4)
    plot(ks4(1,:),p-129,'o-','Color','#a6cee3','DisplayName','Normal','LineWidth',1.5,'MarkerSize',5);
    hold on
    plot(ks4(2,:),p-129,'+--','Color','#1f78b4','DisplayName','Lognormal','LineWidth',1.5,'MarkerSize',5);
    plot(ks4(3,:),p-129,'x-','Color','#b2df8a','DisplayName','Weibull','LineWidth',1.5,'MarkerSize',5);
    plot(ks4(4,:),p-129,'.--','Color','#33a02c','DisplayName','Gamma','LineWidth',1.5,'MarkerSize',5);
    hold off
    set(gca,'YDir','reverse');
    xlabel('p-value');
    ylim([-125 125]);
    title('Lagrangian');

    % Make legend transparent
    h_leg.BoxFace.ColorType='truecoloralpha';
    h_leg.BoxFace.ColorData=uint8(255*[1 1 1 0.75]');

else 
    % WINTER
    subplot(1,4,1)
    plot(ks1(1,:),p,'o-','Color','#a6cee3','DisplayName','Normal','LineWidth',1.5,'MarkerSize',5);
    hold on
    plot(ks1(2,:),p,'+--','Color','#1f78b4','DisplayName','Lognormal','LineWidth',1.5,'MarkerSize',5);
    plot(ks1(3,:),p,'x-','Color','#b2df8a','DisplayName','Weibull','LineWidth',1.5,'MarkerSize',5);
    plot(ks1(4,:),p,'.--','Color','#33a02c','DisplayName','Gamma','LineWidth',1.5,'MarkerSize',5);
    if nargin > 6
        yline(MLD(1),'Color',[0.4 0.4 0.4],'LineWidth',3,'DisplayName','Average MLD','HandleVisibility','off');
    end
    if nargin > 7
        yline(dcmMu(1),'Color',[0.4 0.4 0.4],'LineWidth',3,'DisplayName','DCM','HandleVisibility','off');
    end
    hold off
    set(gca,'YDir','reverse');
    ylim([min(p) max(p)]);
    xlabel('p-value');
    ylabel('Pressure [db]');
    title('Winter');
    
    % SPRING
    subplot(1,4,2)
    plot(ks2(1,:),p,'o-','Color','#a6cee3','DisplayName','Normal','LineWidth',1.5,'MarkerSize',5);
    hold on
    plot(ks2(2,:),p,'+--','Color','#1f78b4','DisplayName','Lognormal','LineWidth',1.5,'MarkerSize',5);
    plot(ks2(3,:),p,'x-','Color','#b2df8a','DisplayName','Weibull','LineWidth',1.5,'MarkerSize',5);
    plot(ks2(4,:),p,'.--','Color','#33a02c','DisplayName','Gamma','LineWidth',1.5,'MarkerSize',5);
    if nargin > 6
        yline(MLD(2),'Color',[0.4 0.4 0.4],'LineWidth',3,'DisplayName','Average MLD','HandleVisibility','off');
    end
    if nargin > 7
        yline(dcmMu(2),'Color',[0.4 0.4 0.4],'LineWidth',3,'DisplayName','DCM','HandleVisibility','off');
    end
    hold off
    set(gca,'YDir','reverse');
    ylim([min(p) max(p)]);
    xlabel('p-value');
    title('Spring');
    
    % SUMMER
    subplot(1,4,3)
    plot(ks3(1,:),p,'o-','Color','#a6cee3','DisplayName','Normal','LineWidth',1.5,'MarkerSize',5);
    hold on
    plot(ks3(2,:),p,'+--','Color','#1f78b4','DisplayName','Lognormal','LineWidth',1.5,'MarkerSize',5);
    plot(ks3(3,:),p,'x-','Color','#b2df8a','DisplayName','Weibull','LineWidth',1.5,'MarkerSize',5);
    plot(ks3(4,:),p,'.--','Color','#33a02c','DisplayName','Gamma','LineWidth',1.5,'MarkerSize',5);
    if nargin > 6
        yline(MLD(3),'Color',[0.4 0.4 0.4],'LineWidth',3,'DisplayName','Average MLD','HandleVisibility','off');
    end
    if nargin > 7
        yline(dcmMu(3),'Color',[0.4 0.4 0.4],'LineWidth',3,'DisplayName','DCM','HandleVisibility','off');
    end
    hold off
    h_leg = legend('Location','best');
    set(gca,'YDir','reverse');
    xlabel('p-value');
    ylim([min(p) max(p)]);
    title('Summer');
    
    % AUTUMN
    subplot(1,4,4)
    plot(ks4(1,:),p,'o-','Color','#a6cee3','DisplayName','Normal','LineWidth',1.5,'MarkerSize',5);
    hold on
    plot(ks4(2,:),p,'+--','Color','#1f78b4','DisplayName','Lognormal','LineWidth',1.5,'MarkerSize',5);
    plot(ks4(3,:),p,'x-','Color','#b2df8a','DisplayName','Weibull','LineWidth',1.5,'MarkerSize',5);
    plot(ks4(4,:),p,'.--','Color','#33a02c','DisplayName','Gamma','LineWidth',1.5,'MarkerSize',5);
    if nargin > 6
        yline(MLD(4),'Color',[0.4 0.4 0.4],'LineWidth',3,'DisplayName','Average MLD','HandleVisibility','off');
    end
    if nargin > 7
        yline(dcmMu(4),'Color',[0.4 0.4 0.4],'LineWidth',3,'DisplayName','DCM','HandleVisibility','off');
    end
    hold off
    set(gca,'YDir','reverse');
    xlabel('p-value');
    ylim([min(p) max(p)]);
    title('Autumn');
    
    % Make legend transparent
    h_leg.BoxFace.ColorType='truecoloralpha';
    h_leg.BoxFace.ColorData=uint8(255*[1 1 1 0.75]');
end

end