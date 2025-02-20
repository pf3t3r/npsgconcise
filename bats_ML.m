clear;clc;close all;
format longG;

% omnibus script to recreate all Stn ALOHA figures with BATS data.

showOtherTests = false;     % Evaluate data with Lilliefors and chi^2 tests.
                            % (in addition to Anderson-Darling)
showL0title = false;        % Leave as 'false' for paper.
%% BATS. chl-a. Hovmoeller. (Fig 1a)

% Import Data.
tmp = importdata('data\L0\bats_pigments.txt').data;
id = tmp(:,1);      % bottle ID with format !####$$$@@, where ! = cruise
                    % type (1 = core cruise), #### = cruise number, 
                    % $$$ = cast number, @@ = Niskin number.
YMD = tmp(:,2);     % Year, Month, and Day.
hhmm = tmp(:,4);    % hours and minutes.
QF = tmp(:,7);      % Quality control factor.
depth = tmp(:,8);           % Depth (m).
chla = tmp(:,22);           % chl-a concentration (ng/l).
chla_Turner = tmp(:,24);    % as above but using Turner method (ng/l).

% Re-assign NaNs.
chla(chla==-999) = nan;
chla_Turner(chla_Turner==-999) = nan;

% Set up time vector.
for i = 1:length(YMD)
    tmp1 = num2str(YMD(i));
    YY(i) = str2num(tmp1(1:4));
    MM(i) = str2num(tmp1(5:6));
    DD(i) = str2num(tmp1(7:8));
    clear tmp1;

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
clear YY MM DD hh mm ss;

% Bin by depth.
edges = 0:10:200;
depth_B = discretize(depth,edges);
% depth_B(isnan(depth_B)) = 99;
clear edges;

% Extract cruise type 'cType' and cruise number 'CRN'.
tmp2 = num2str(id);
cType = str2num(tmp2(:,1));
CRN = str2num(tmp2(:,2:5));
cast = str2num(tmp2(:,6:8));
clear tmp2;

% Extract core cruises. These are labelled as cType = 1.
ids = find(cType==1);
CoreCRN = CRN(ids);
cast = cast(ids);
copyOfCoreCrn = CoreCRN;
t = t(ids);
chla = chla(ids);
depth = depth(ids);
depth_B = depth_B(ids);

% Save ID where new cruise starts.
id_nc = [1];
CRN_no = [1];
for i = 2:length(CoreCRN)
    if CoreCRN(i) > CoreCRN(i-1)
        disp(i)
        id_nc = [id_nc i];
        CRN_no = [CRN_no CoreCRN(i)];
    end
end

crnAndCast = [CoreCRN cast];

% Divide up the depth, chlorophyll, and time arrays by cruise.
% NOTE that there are 405 cruises but only 398 have chl-a data.
dgrid = nan(398,35);

% For cruise #1.
dgrid(1,1:(id_nc(2)-id_nc(1))) = depth_B(id_nc(1):id_nc(2)-1);
depthUnbinned(1,1:(id_nc(2)-id_nc(1))) = depth(id_nc(1):id_nc(2)-1);
tgrid(1,1:(id_nc(2)-id_nc(1))) = t(id_nc(1):id_nc(2)-1);
chlagrid(1,1:(id_nc(2)-id_nc(1))) = chla(id_nc(1):id_nc(2)-1);

% Loop for cruises #2-397.
for i = 2:397
    dgrid(i,1:(id_nc(i+1)-id_nc(i))) = depth_B(id_nc(i):id_nc(i+1)-1);
    depthUnbinned(i,1:(id_nc(i+1)-id_nc(i))) = depth(id_nc(i):id_nc(i+1)-1);
    tgrid(i,1:(id_nc(i+1)-id_nc(i))) = t(id_nc(i):id_nc(i+1)-1);
    chlagrid(i,1:(id_nc(i+1)-id_nc(i))) = chla(id_nc(i):id_nc(i+1)-1);
    %pgrid(3,1:(newCruiseId(4)-newCruiseId(3))) = pB(newCruiseId(3):newCruiseId(4)-1);
end

% Final cruise (#398).
dgrid(398,1:(length(depth_B)+1-id_nc(398))) = depth_B(id_nc(398):length(depth_B));
depthUnbinned(398,1:(length(depth)+1-id_nc(398))) = depth(id_nc(398):length(depth));
tgrid(398,1:(length(t)+1-id_nc(398))) = t(id_nc(398):length(t));
chlagrid(398,1:(length(chla)+1-id_nc(398))) = chla(id_nc(398):length(chla));

tgridDatenum = datenum(tgrid);      % Convert time vector to datenum format (for plotting).

% Figure 1a. chl-a(p,t).
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

%% BATS. chla-a. Mean profile. (Fig. 1b)

% Calculate the mean chl-a at each depth bin as well as the 5th and 95th
% percentile values.
chlaProfile = nan(1,20); f5 = nan(1,20); f95 = nan(1,20);
for i = 1:20
    chlaProfile(i) = mean(chlagrid(dgrid==i),"omitnan");
    f5(i) = prctile(chlagrid(dgrid==i),5);
    f95(i) = prctile(chlagrid(dgrid==i),95);
end

% Toggle show y-label and title (for the paper we don't need either)
displayYLabelAndTitle = false;

% Figure 1b. Mean chl-a profile.
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
ylim([1 18]); xlim([0 600]);
ax = gca;
ax.FontSize = 15;


%% BATS. chl-a. A-D test. (Fig. 2)

% Find the depth of the Chlorophyll Maximum (CM). NOTE that in the case of
% BATS that this is not necessarily a deep maximum. We will separate the
% deep maxima at a later stage.
depthOfCm = nan(398,1);
for i = 1:398
    tmp = chlagrid(i,:);
    [~,id_DCM(i)] = max(tmp);
    depthOfCm(i) = depthUnbinned(i,id_DCM(i));
end

% Calculate the mean depth of the Chlorophyll Maximum (CM) as well as the
% 5th and 9th percentile interval values.
meanCM = mean(depthOfCm,"omitnan");
CM_5pct = prctile(depthOfCm,5);
CM_95pct = prctile(depthOfCm,95);

% Use the Anderson-Darling (A-D) test to evaluate whether the data is
% distributed normally or lognormally. The hypothesis test result 'h' will
% return as h = 1 if the null hypothesis is rejected or h = 0 if there is a
% failure to reject the null hypothesis.
hN = nan(20,1); pN = nan(20,1); hL = nan(20,1); pL = nan(20,1);
obs = nan(20,1);
for i = 1:20
    tmp = chla(depth_B==i);
    if length(tmp) > 30
        obs(i) = length(tmp);
        tmp(tmp==0) = nan;

        [hN(i), pN(i)] = adtest(tmp,"Distribution","norm","Alpha",0.005);
        [hL(i), pL(i)] = adtest(tmp,"Distribution","logn","Alpha",0.005);

        if showOtherTests == true
            pd0 = fitdist(tmp,'Normal');
            pd = fitdist(tmp,'Lognormal');
            [hN2(i), pN2(i)] = chi2gof(tmp,"CDF",pd0);
            [hN3(i), pN3(i)] = lillietest(tmp,"Distribution","norm");
            [hx1(i),px1(i)] = chi2gof(tmp,'CDF',pd);
            [hx2(i),px2(i)] = lillietest(log(tmp),"Distr","norm");
        end
    end
end

% Figure 2. chl-a. Is the data normal or lognormal?
figure;
if showL0title == true
    sgtitle("chl-a (L0): " + "BATS "+num2str(YMD(1))+" - " + num2str(YMD(end))+"");
end
subplot(1,3,[1 2])
yyaxis left
semilogx(pN,0.5:1:19.5,'o-','Color','#c51b7d','DisplayName','Normal (A-D)','LineWidth',1.5,'MarkerSize',5);
hold on
semilogx(pL,0.5:1:19.5,'o-','Color','#4d9221','DisplayName','Lognormal (A-D)','LineWidth',1.5,'MarkerSize',5);
if showOtherTests == true
    semilogx(pN2,1:1:20,'o-','Color','#c51b7d','DisplayName','Normal (chi^2)','LineWidth',1.5,'MarkerSize',5);
    semilogx(pN3,1:1:20,'o--','Color','#c51b7d','DisplayName','Normal (Lil.)','LineWidth',1.5,'MarkerSize',5);
    semilogx(px1,1:1:20,'o-','Color','#4d9221','DisplayName','Lognormal (chi^2)','LineWidth',1.5,'MarkerSize',5);
    semilogx(px2,1:1:20,'o--','Color','#4d9221','DisplayName','Lognormal (Lil.)','LineWidth',1.5,'MarkerSize',5);
end
set(gca,"YDir","reverse"); legend();
yticklabels(0:20:200);
ylim([0 20]);
ylabel("Depth (m) (10-m bins)");
yyaxis right
yline(meanCM,DisplayName="p_{DCM} \pm 5/95",Interpreter="latex");
yline(CM_5pct,HandleVisibility="off");
yline(CM_95pct,HandleVisibility="off");
xline(0.005,":","\alpha","DisplayName","\alpha = 0.005");
hold off
set(gca,"YDir","reverse"); legend();
yticklabels({});
ylim([0 200]);
xlim([1e-3 1]);
xlabel("p-value");

subplot(1,3,3)
barh(obs,'FaceColor','#d3d3d3');
hold on
xline(30);
set(gca,"YDir","reverse"); xlabel("No. of Obs.");
ylim([0.5 20.5]); yticklabels({});

%% Calculate MLD for BATS from temperature and salinity data. (prep Fig. 3)

% Import data.
data = importdata("data\bats_bottle.txt").data;
id = data(:,1);     % Bottle ID with format !####$$$@@, where ! = cruise
                    % type (1 = core cruise), #### = cruise number, 
                    % $$$ = cast number, @@ = Niskin number.
YMD = data(:,2);                    % Year, Month, Day.
lat = data(:,5); long = data(:,6);  % Latitude and Longitude.
z = data(:,8);                      % Depth (m).
T = data(:,9);                      % Temperature (C).
Sp = data(:,11);                    % Practical Salinity.

% Setup time vector.
tmpB = num2str(YMD);
YY = str2num(tmpB(:,1:4));
MM = str2num(tmpB(:,5:6));
DD = str2num(tmpB(:,7:8));
t = datetime(YY,MM,DD);

% (re-) assign NaNs.
z(z==-999) = nan;

% Process Hydrography data according to TEOS-10 standard.
p = gsw_p_from_z(-z,lat);
SA = gsw_SA_from_SP(Sp,p,long,lat);
CT = gsw_CT_from_t(SA,T,p);

% Import cruise type "cType" and cruise number "CRN". 
tmp3 = num2str(id);
cType = str2num(tmp3(:,1));
CRN = str2num(tmp3(:,2:5));
clear tmp3;

% Extract hydrography from the core cruises (labelled as cType = 1).
ids = find(cType==1);
CoreCRN = CRN(ids);
SA = SA(ids);
CT = CT(ids);
p = p(ids);

% Save the point at which a new cruise starts to the array 'id_nc'.
% Additionally, save the actual cruise number of the new cruise to the
% array 'crn_mr'. This is necessary because not all cruises have MLD data.
id_nc = [1];    % ID (id_) of new cruise (nc).
crn_mr = [1];   % cruise number (crn_) where MLD was recorded (mr).
for i = 2:length(CoreCRN)
    if CoreCRN(i) > CoreCRN(i-1)
        id_nc = [id_nc i];
        crn_mr = [crn_mr CoreCRN(i)];
    end
end
t_nc = t(id_nc);    % Time at which new cruise starts.

% Calculate the Mixed Layer Depth per cruise 'MLD_pc'
mld_pc = 20*ones(404,1);        % Minimum possible MLD = 20 dbar.
for i = 2:400
    mld_pc(i) = gsw_mlp(SA(id_nc(i):id_nc(i+1)),CT(id_nc(i):id_nc(i+1)),p(id_nc(i):id_nc(i+1)));
end

% Remove casts that are unrealistic.
mld_pc(mld_pc>300) = nan;
mld_pc(isnan(mld_pc)) = 20;

%% DCM vs MLP. In terms of CRN.

% Remove unrealistic values.
id_DCM(id_DCM>300)=nan;

% Convert DCM depth (m) to DCM pressure (dbar).
p_DCM = gsw_p_from_z(-depthOfCm,lat(1));

figure;
plot(CRN_no,p_DCM,DisplayName="DCM");
hold on
plot(crn_mr,mld_pc,DisplayName="MLD");
hold off
set(gca,"YDir","reverse");
ylabel("Pressure [dbar]"); xlabel("Cruise No.");
legend();


% Which cruises to look at with Stn ALOHA Level Methodology. We need the
% DCM to be beneath the ML.
% Consolidate two vectors of MLD and DCM so they have the same size.

cruises = 1:1:405;                          % number per cruise
newDCMcrnVector = nan(length(cruises),1);
newMLDcrnVector = nan(length(cruises),1);
newMlp = nan(length(cruises),1);
newDcm = nan(length(cruises),1);

for i = 1:405
    for j = 1:398
        if CRN_no(j) == i
            newDCMcrnVector(i) = CRN_no(j);
            newDcm(i) = p_DCM(j);
        end
    end
    for k = 1:404
        if crn_mr(k) == i
            newMLDcrnVector(i) = crn_mr(k);
            newMlp(i) = mld_pc(k);
        end
    end
end

figure;
plot(newDCMcrnVector,newDcm,DisplayName="Chl-a Maximum");
hold on
plot(newMLDcrnVector,newMlp,DisplayName="Mixed Layer Depth");
hold off
set(gca,"YDir","reverse");
ylabel("Pressure [dbar]"); xlabel("Cruise No.");
legend();

%% Check where DCM is beneath MLD.

dcmBelow = newDcm - newMlp;
cruisesWhereDCMisBelowMLD = [];

for i = 1:405
    if dcmBelow(i) > 0
        cruisesWhereDCMisBelowMLD = [cruisesWhereDCMisBelowMLD i]
    end
end


figure
plot(1:1:405,dcmBelow);

% add column to start
newChla = [CRN_no' chlagrid];

test1 = 1:1:405;
test1 = test1';
tmpC = nan(length(test1),35); tmpD = nan(length(test1),35); tmpT = nan(length(test1),35);
for i = 1:405
    for j = 1:398
        if CRN_no(j) == i
            tmpC(i,:) = chlagrid(j,:);
            tmpD(i,:) = dgrid(j,:);
            tmpT(i,:) = datenum(tgrid(j,:));
        end
    end
end

test3 = [test1 tmpC]; % this has forced the chlagrid array onto a 405 cruise matrix.
test4 = [test1 tmpD]; % same for depth
test5 = [test1 tmpT]; % ... and time
% Now show only the chla arrays where DCM is beneath the MLD.
newChlaArray = test3(cruisesWhereDCMisBelowMLD,:);
newDepthArray = test4(cruisesWhereDCMisBelowMLD,:);
newTimeArray = test5(cruisesWhereDCMisBelowMLD,:);

% remember that the first entry is the CRUISE NO.

dcmsBelow = newDcm(cruisesWhereDCMisBelowMLD);

% try to work out which bottles are from these cruises where the DCM is
% below the MLD.

test66 = [];
for i = 1:length(copyOfCoreCrn)
    for j = 1:length(cruisesWhereDCMisBelowMLD)
        if copyOfCoreCrn(i) == cruisesWhereDCMisBelowMLD(j)
            test66 = [test66 i];
        end
    end
end

chla_lowDCM = chla(test66);
dep_lowDCM = depth_B(test66);

depthOfCm = nan(267,1);
for i = 1:267
    tmp = newChlaArray(i,:);
    [~,id_DCM(i)] = max(tmp);
    %depthOfCm(i) = depthUnbinned(i,id_DCM(i));
    depthOfCm(i) = id_DCM(i);
end

depthOfCm = 10*depthOfCm - 5;

% Calculate the mean depth of the Chlorophyll Maximum (CM) as well as the
% 5th and 9th percentile interval values.
meanCM = mean(depthOfCm,"omitnan");
CM_5pct = prctile(depthOfCm,5);
CM_95pct = prctile(depthOfCm,95);


% Use the Anderson-Darling (A-D) test to evaluate whether the data is
% distributed normally or lognormally. The hypothesis test result 'h' will
% return as h = 1 if the null hypothesis is rejected or h = 0 if there is a
% failure to reject the null hypothesis.
hN = nan(20,1); pN = nan(20,1); hL = nan(20,1); pL = nan(20,1);
obs = nan(20,1);
for i = 1:20
    tmp = chla_lowDCM(dep_lowDCM==i);
    if length(tmp) > 30
        obs(i) = length(tmp);
        tmp(tmp==0) = nan;

        [hN(i), pN(i)] = adtest(tmp,"Distribution","norm","Alpha",0.005);
        [hL(i), pL(i)] = adtest(tmp,"Distribution","logn","Alpha",0.005);

        if showOtherTests == true
            pd0 = fitdist(tmp,'Normal');
            pd = fitdist(tmp,'Lognormal');
            [hN2(i), pN2(i)] = chi2gof(tmp,"CDF",pd0);
            [hN3(i), pN3(i)] = lillietest(tmp,"Distribution","norm");
            [hx1(i),px1(i)] = chi2gof(tmp,'CDF',pd);
            [hx2(i),px2(i)] = lillietest(log(tmp),"Distr","norm");
        end
    end
end

% Figure 2X. chl-a. L0. Is the data normal or lognormal? This time we look
% at only those cruises where the DCM was beneath the MLD.

figure;
if showL0title == true
    sgtitle("chl-a (L0): " + "BATS "+num2str(YMD(1))+" - " + num2str(YMD(end))+"");
end
subplot(1,3,[1 2])
yyaxis left
semilogx(pN,0.5:1:19.5,'o-','Color','#c51b7d','DisplayName','Normal (A-D)','LineWidth',1.5,'MarkerSize',5);
hold on
semilogx(pL,0.5:1:19.5,'o-','Color','#4d9221','DisplayName','Lognormal (A-D)','LineWidth',1.5,'MarkerSize',5);
if showOtherTests == true
    semilogx(pN2,1:1:20,'o-','Color','#c51b7d','DisplayName','Normal (chi^2)','LineWidth',1.5,'MarkerSize',5);
    semilogx(pN3,1:1:20,'o--','Color','#c51b7d','DisplayName','Normal (Lil.)','LineWidth',1.5,'MarkerSize',5);
    semilogx(px1,1:1:20,'o-','Color','#4d9221','DisplayName','Lognormal (chi^2)','LineWidth',1.5,'MarkerSize',5);
    semilogx(px2,1:1:20,'o--','Color','#4d9221','DisplayName','Lognormal (Lil.)','LineWidth',1.5,'MarkerSize',5);
end
set(gca,"YDir","reverse"); legend();
yticklabels(0:20:200);
ylim([0 20]);
ylabel("Depth (m) (10-m bins)");
yyaxis right
yline(meanCM,DisplayName="p_{DCM} \pm 5/95",Interpreter="latex");
yline(CM_5pct,HandleVisibility="off");
yline(CM_95pct,HandleVisibility="off");
xline(0.005,":","\alpha","DisplayName","\alpha = 0.005");
hold off
set(gca,"YDir","reverse"); legend();
yticklabels({});
ylim([0 200]);
xlim([1e-3 1]);
xlabel("p-value");

subplot(1,3,3)
barh(obs,'FaceColor','#d3d3d3');
hold on
xline(30);
set(gca,"YDir","reverse"); xlabel("No. of Obs.");
ylim([0.5 20.5]); yticklabels({});

%% Find DCM according to CTD.

