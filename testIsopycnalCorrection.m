close all;clc;clear;

ctd = load("datafiles\ctd_iso_ALL.mat").ctd;
ct50_ctd = ctd(50).ct;
t50 = datetime(ctd(50).date,"ConvertFrom","datenum");
p = ctd(50).p;


%%
% iso = load("datafiles\ctd_iso_ALL.mat").iso;
% ct50_iso = iso(50).ct;
% deltaT = ct50_ctd - ct50_iso;

% plotThis = false;
% 
% if plotThis
%     figure;
%     subplot(1,2,1)
%     plot(ct50_ctd(:,5),p,DisplayName="\Theta");
%     hold on
%     plot(ct50_iso(:,5),p,DisplayName="\Theta_{iso}");
%     hold off
%     set(gca,"YDir","reverse");
%     title("$\Theta$ vs. $\Theta_{iso}$","Interpreter","latex")
%     sgtitle("HOT-50: Cast 5","FontSize",8);
%     legend(Location="southeast"); grid on;
%     ylabel("Pressure [dbar]");
%     xlabel("\Theta [{\circ}C]");
%     
%     subplot(1,2,2)
%     plot(deltaT(:,5),p);
%     set(gca,"YDir","reverse");
%     title("$\Delta \Theta = \Theta - \Theta_{iso}$","Interpreter","latex")
%     yticklabels({}); grid on;
%     xlabel("\Theta [{\circ}C]");
% end

%%
% clear deltaT; 
sig50_ctd = ctd(50).sig(1:101,:);
p = p(1:101,1);

% set up binned densities
sbm = floor(100*min(min(sig50_ctd)))./100;
sbM = ceil(100*max(max(sig50_ctd)))./100;
step = 0.05;

% binned densities
sb = sbm:step:sbM;

n = length(ct50_ctd(1,:));
meanPressurePerIsopycnal = nan(length(sb),n+2);
% 0. Loop for each cast
for l = 1:n
    disp(n);

    % 1. find pressures that correspond to each binned density
    depthBins = nan(101,3);
    for i = 1:101
        for j = 1:length(sb)
            if j < length(sb)
                if sig50_ctd(i,l) > sb(j) & sig50_ctd(i,l) < sb(j+1)
                    depthBins(i,1) = sb(j);
                    depthBins(i,2) = j;
                    depthBins(i,3) = p(i);
                end
            else
                if sig50_ctd(i,l) > sb(j)
                    depthBins(i,1) = sb(j);
                    depthBins(i,2) = j;
                    depthBins(i,3) = p(i);
                end
            end
        end
    end
    
    newIsopycnal = [1];
    for i = 2:101
        if depthBins(i,2) > depthBins(i-1,2)
            newIsopycnal = [newIsopycnal i];
        end
    end
    
    distinctBins = unique(depthBins(:,2));
    
    meanPressure = [];
    for k = 2:length(newIsopycnal)
        meanPressure = [meanPressure mean(depthBins(newIsopycnal(k-1):newIsopycnal(k)-1,3))];
    end
    meanPressure = [meanPressure mean(depthBins(newIsopycnal(end):end,3))];
    
    meanPressurePerBin = [distinctBins meanPressure'];
    
    % 2. find mean pressure at that density for that cast
    
    for k = 1:length(sb)
        meanPressurePerIsopycnal(k,1) = sb(k);
        meanPressurePerIsopycnal(k,2) = k;
    end
    
    meanPressurePerIsopycnal(meanPressurePerBin(:,1),l+2) = meanPressurePerBin(:,2);

    clear depthBins meanPressurePerBin i j k l step;
end

figure
plot(meanPressurePerIsopycnal(:,3:end));
xticks(1:1:51); xticklabels(sb);
xlabel("\sigma_0 [kg m^{-3}]"); ylabel("Pressure [dbar]");
set(gca,"YDir","reverse");
legend();

% Depth of sigma = 25 kg/m3 isopycnal (defined as 24.00 - 24.05)
figure
plot(t50,meanPressurePerIsopycnal(30,3:end),'o-',DisplayName="24.60 kgm^{-3}");
hold on
plot(t50,meanPressurePerIsopycnal(34,3:end),'o-',DisplayName="24.80 kgm^{-3}");
plot(t50,meanPressurePerIsopycnal(38,3:end),'o-',DisplayName="25.00 kgm^{-3}");
plot(t50,meanPressurePerIsopycnal(42,3:end),'o-',DisplayName="25.20 kgm^{-3}");
hold off
set(gca,"YDir","reverse");
xlabel("time");
ylabel("pressure [dbar]");
title("Isopycnal Depth: HOT-50");
grid on; grid minor;
legend();

%% Pressure range per isopycnal

for i = 1:51
    pm(i) = min(meanPressurePerIsopycnal(i,3:end));
    pM(i) = max(meanPressurePerIsopycnal(i,3:end));
    pmn(i) = mean(meanPressurePerIsopycnal(i,3:end),"omitnan");
    Dp(i) = pM(i) - pm(i);
end

a = string(round(pmn));
f5s = 58*ones(51,1); % array of 58s: just for plotting text labels

figure;
plot(Dp,sb,'square-','Color',[0.4940 0.1840 0.5560]);
text(f5s,sb,a,"Color",[0.5 0.5 0.5],"FontSize",8);
set(gca,"YDir","reverse");
xlabel('\Delta p = p_{max} - p_{min} [dbar]'); 
ylabel("\sigma_0 [kg m^{-3}]");
title("Pressure Range Across Isopycnals: HOT-50",FontSize=13);
subtitle("28-30/10/1993 (18 casts; p_{mean} in grey)","FontSize",9);
grid on; grid minor;



