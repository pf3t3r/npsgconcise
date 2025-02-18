clear;clc;close all;
format longG;

% omnibus script to recreate all Stn ALOHA figures with BATS data.

%% BATS. chl-a. Hovmoeller. (Fig 1a)

tmp = importdata('data\L0\bats_pigments.txt').data;

YMD = tmp(:,2);
hhmm = tmp(:,4);
for i = 1:length(YMD)
    tmpA = num2str(YMD(i));
    YY(i) = str2num(tmpA(1:4));
    MM(i) = str2num(tmpA(5:6));
    DD(i) = str2num(tmpA(7:8));

    tmpB = num2str(hhmm(i));
    if length(tmpB) == 3
        hh(i) = str2num(tmpB(1));
        mm(i) = str2num(tmpB(2:3));
    elseif length(tmpB) == 4
        hh(i) = str2num(tmpB(1:2));
        mm(i) = str2num(tmpB(3:4));
    else
        hh(i) = nan;
        mm(i) = nan;
    end
end
ss = zeros(1,length(mm));
t = datetime(YY,MM,DD,hh,mm,ss);

QF = tmp(:,7);
depth = tmp(:,8);
chla = tmp(:,22);
chla_Turner = tmp(:,24);
chla(chla==-999) = nan;
chla_Turner(chla_Turner==-999) = nan;

% bin by depth
edges = 0:10:200;
depth_B = discretize(depth,edges);
depth_B(isnan(depth_B)) = 99;

newCast = [1];

for i = 2:length(depth_B)
    if depth_B(i) < depth_B(i-1)
        disp(i);
        newCast = [newCast i];
    end
end

dgrid = nan(647,20);
tgrid = NaT(647,20);
chlagrid = NaN(647,20);

% First Case.
dgrid(1,1:(newCast(2)-newCast(1))) = depth_B(newCast(1):newCast(2)-1);
tgrid(1,1:(newCast(2)-newCast(1))) = t(newCast(1):newCast(2)-1);
chlagrid(1,1:(newCast(2)-newCast(1))) = chla(newCast(1):newCast(2)-1);
% Then Loop.
% i = 2:319
for i = 2:646
    dgrid(i,1:(newCast(i+1)-newCast(i))) = depth_B(newCast(i):newCast(i+1)-1);
    tgrid(i,1:(newCast(i+1)-newCast(i))) = t(newCast(i):newCast(i+1)-1);
    chlagrid(i,1:(newCast(i+1)-newCast(i))) = chla(newCast(i):newCast(i+1)-1);
    %pgrid(3,1:(newCast(4)-newCast(3))) = pB(newCast(3):newCast(4)-1);
end

% Final Case.
dgrid(647,1:(length(depth_B)+1-newCast(647))) = depth_B(newCast(647):length(depth_B));
tgrid(647,1:(length(t)+1-newCast(647))) = t(newCast(647):length(t));
chlagrid(647,1:(length(chla)+1-newCast(647))) = chla(newCast(647):length(chla));

tgridDatenum = datenum(tgrid);

figure;
contourf(tgridDatenum,dgrid,chlagrid,linspace(0,500,10),'LineColor','auto');
set(gca,"YDir","reverse");
datetickzoom('x','yyyy','keeplimits');
colormap(flipud(cbrewer2("GnBu")));
c = colorbar;
c.Label.String = 'chl-a [ng/l]';
c.FontSize = 13;
ylim([1 18]);
zlim([0 500]);
yticks(1:1:18);
% yticklabels({});
yticklabels(5:10:175);
ylabel("P [dbar]","FontSize",13); xlabel("Time",FontSize=13);
ax = gca;
ax.FontSize = 15;

% Filter by season.
tgridF = tgrid;
chlagridF = chlagrid;
for i = 1:647
    for j = 1:20
        if month(tgridF(i,j)) == 12 || month(tgridF(i,j)) == 1 || month(tgridF(i,j)) == 2
            chlagridF(i,j) = NaN;
        else
            chlagridF(i,j) = chlagrid(i,j);
        end
    end
end

tgridFD = datenum(tgridF);

figure;
contourf(tgridFD,dgrid,chlagridF,linspace(0,500,10),'LineColor','auto');
set(gca,"YDir","reverse");
datetickzoom('x','yyyy','keeplimits');
colormap(flipud(cbrewer2("GnBu")));
c = colorbar;
c.Label.String = 'chl-a [ng/l]';
c.FontSize = 13;
ylim([1 18]);
zlim([0 500]);
yticks(1:1:18);
% yticklabels({});
yticklabels(5:10:175);
ylabel("P [dbar]","FontSize",13); xlabel("Time",FontSize=13);
ax = gca;
ax.FontSize = 15;

%% BATS. chla-a. Mean profile. (Fig. 1b)
chlaProfile = nan(1,20);
f5 = nan(1,20);
f95 = nan(1,20);

for i = 1:20
    chlaProfile(i) = mean(chlagrid(dgrid==i),"omitnan");
    f5(i) = prctile(chlagrid(dgrid==i),5);
    f95(i) = prctile(chlagrid(dgrid==i),95);
end

% Toggle show y-label and title (for paper we don't need either)
displayYLabelAndTitle = false;

% Plot the mean profile of fluorescence with the mean and confidence
% interval.
ax = figure;
plot(chlagrid(:,1),dgrid(:,1),'.',Color=[0.8 0.8 0.8],DisplayName="raw data");
hold on
plot(chlagrid(:,2:20),dgrid(:,2:20),'.',Color=[0.8 0.8 0.8],HandleVisibility='off');
plot(chlaProfile,1:1:20,'-',"Color",[0 0 0],DisplayName="mean");
plot(f5,1:1:20,'-',"Color",[0.5 0.5 0.5],DisplayName="5%");
plot(f95,1:1:20,'-',"Color",[0.5 0.5 0.5],DisplayName="95%");
hold off
set(gca,"YDir","reverse");
legend();
xlabel("chl-$a$ [ng/L]",Interpreter="latex");
if displayYLabelAndTitle == true
    title("L0 Chl-$a$ 1988-2022",Interpreter="latex");
    ylabel("P [dbar]",Interpreter="latex");
    yticklabels({});
end
yticks(1:1:18); yticklabels(5:10:175);
ylim([1 18]);
ax = gca;
ax.FontSize = 15;


%% BATS. chl-a. A-D test. (Fig. 2)

% mean DCM and prctl.
for i = 1:647
    tmp = chlagridF(i,:);
    [mx(i) pcm(i)] = max(tmp);
    if pcm(i) == 1 || pcm(i) == 2
        pcm(i) = nan;
        mx(i) = nan;
    end
end

% ppp =id*10 -5;
ppp2 = mean(pcm,"omitnan");
prct1 = prctile(pcm,5);
prct2 = prctile(pcm,95);

fluo_B = nan(20,1000);
hN = nan(20,1); pN = nan(20,1); hL = nan(20,1); pL = nan(20,1);
obs = nan(20,1);
for i = 1:20
    tmp = chla(depth_B==i);
    if length(tmp) > 30
        obs(i) = length(tmp);
        tmp(tmp==0) = nan;

        [hN(i), pN(i)] = adtest(tmp,"Distribution","norm");
        [hL(i), pL(i)] = adtest(tmp,"Distribution","logn","Alpha",0.005);

        
        pd0 = fitdist(tmp,'Normal');
        pd = fitdist(tmp,'Lognormal');
        [hN2(i), pN2(i)] = chi2gof(tmp,"CDF",pd0);
        [hN3(i), pN3(i)] = lillietest(tmp,"Distribution","norm");
        [hx1(i),px1(i)] = chi2gof(tmp,'CDF',pd);
        [hx2(i),px2(i)] = lillietest(log(tmp),"Distr","norm");
    end
end

figure;
sgtitle("chl-a (L0): " + "BATS "+num2str(YMD(1))+" - " + num2str(YMD(end))+"");
subplot(1,2,1)
semilogx(pN,1:1:20,'o-','Color','#c51b7d','DisplayName','Normal (A-D)','LineWidth',1.5,'MarkerSize',5);
hold on
% semilogx(pN2,1:1:20,'o-','Color','#c51b7d','DisplayName','Normal (chi^2)','LineWidth',1.5,'MarkerSize',5);
% semilogx(pN3,1:1:20,'o--','Color','#c51b7d','DisplayName','Normal (Lil.)','LineWidth',1.5,'MarkerSize',5);
semilogx(pL,1:1:20,'o-','Color','#4d9221','DisplayName','Lognormal (A-D)','LineWidth',1.5,'MarkerSize',5);
% semilogx(px1,1:1:20,'o-','Color','#4d9221','DisplayName','Lognormal (chi^2)','LineWidth',1.5,'MarkerSize',5);
% semilogx(px2,1:1:20,'o--','Color','#4d9221','DisplayName','Lognormal (Lil.)','LineWidth',1.5,'MarkerSize',5);
yline(ppp2,DisplayName="p_{DCM} \pm 5/95",Interpreter="latex");
yline(prct1,HandleVisibility="off");
yline(prct2,HandleVisibility="off");
xline(0.005,":","\alpha","DisplayName","\alpha = 0.005");
hold off
set(gca,"YDir","reverse"); legend();
yticklabels(-5:20:195);
ylim([0.5 20.5]);
xlim([1e-3 1]);
ylabel("Depth (m) (10-m bins)");
xlabel("p-value");
subplot(1,2,2)
barh(obs,'FaceColor','#d3d3d3');
hold on
xline(30);
set(gca,"YDir","reverse");
ylim([0.5 20.5]); yticklabels({});
title("No. of Obs.");


% show DCM per cruise
figure;
plot(mean(tgridF,2)',pcm); 
set(gca,"YDir","reverse");
yticklabels(-5:20:195);
%% Calculate MLD for BATS from temperature and salinity data. (prep Fig. 3)

data = importdata("data\bats_bottle.txt").data;

id = data(:,1);
lat = data(:,5);
long = data(:,6);
z = data(:,8);
T = data(:,9);
Sp = data(:,11);

z(z==-999) = nan;

% Import Hydrography data.
p = gsw_p_from_z(-z,lat);
SA = gsw_SA_from_SP(Sp,p,long,lat);
CT = gsw_CT_from_t(SA,T,p);


% % Test: find mixed layer pressure for cruise 2.
% a=25; b=72;
% mlp = gsw_mlp(SA(a:b),CT(a:b),p(a:b));

% Try to generalise the above. Import cruise type "cType" and cruise number
% "CRN". 
tmpA = num2str(id);
cType = str2num(tmpA(:,1));
CRN = str2num(tmpA(:,2:5));

% Extract the core cruises (labelled as cType = 1).
ids = find(cType==1);
CoreCRN = CRN(ids);
SA = SA(ids);
CT = CT(ids);
p = p(ids);

% % % Remove cruises 22, 23, and 32. (they give unrealistic MLP)
% CoreCRN = CoreCRN(find(CoreCRN~=22 & CoreCRN~=23 & CoreCRN~=32));
% SA = SA(find(CoreCRN~=22 & CoreCRN~=23 & CoreCRN~=32));
% CT = CT(find(CoreCRN~=22 & CoreCRN~=23 & CoreCRN~=32));
% p = p(find(CoreCRN~=22 & CoreCRN~=23 & CoreCRN~=32));

% Save ID of when new cruise starts.
newCruiseId = [1];
for i = 2:length(CoreCRN)
    if CoreCRN(i) > CoreCRN(i-1)
        disp(i)
        newCruiseId = [newCruiseId i];
    end
end

% find MLP per cruise
mlpByCRN = nan(404,1);
for i = 2:400
    mlpByCRN(i) = gsw_mlp(SA(newCruiseId(i):newCruiseId(i+1)),CT(newCruiseId(i):newCruiseId(i+1)),p(newCruiseId(i):newCruiseId(i+1)));
end
mlpByCRN(mlpByCRN>300) = nan;

figure
plot(mlpByCRN);
