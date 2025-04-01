% Statistical Analysis of upper 200 dbar at Station ALOHA for chl-a, other pigments, and BGC variables.
% We import the data and run a hypothesis test on it with the respective
% null hypotheses of normal and lognormal. We use the Anderson-Darling test
% since this is both more powerful than similar tests such as
% Kolmogorov-Smirnov and more flexible than tests such as Shapiro-Wilks
% which did not easily allow for testing of other distributions.

clear; clc; close all;
addpath("func\");
set(groot, 'defaultFigureUnits', 'centimeters', 'defaultFigurePosition', [3 3 15 15]);

% Possible test cases.
principle = true;       % main analysis
seasonal = false;       % seasonality of statistics: A-D
startYear = false;      % analyse effect (if any) of varying start year on 
                        % distributions
newCtd = false;         % analyse 2001-2021 data (to mirror CTD results)
night = false;          % analyse night-time 2001-2021 (to mirror CTD results)
logAxes = true;         % output p-values as log values
showL0title = true;     % show title for chl-a plot. (off for paper)

if logAxes == true
    lp = "log/";
else
    lp = "";
end

%% Chl-a depth and time-series.
% Import Chl-a data and transform the data.
% Plot the depth and time-series of chla across the upper 200 dbar (approx
% 200 m) between 1988 and 2021.

% Import data file
tmp = importdata('data/L0/hplcChla_88-21_200.txt');

% Import binned pressure and concentration of chl-a
[~,~,~,pB,chla,~] = L0_helper(tmp,50,'ad',true,0,true);

time = tmp.data(:,2);
hms = tmp.data(:,3);

YY1 = mod(time,100);

% Correct the years
for i=1:length(YY1)
    if YY1(i) < 85
        YY1(i) = YY1(i) + 2000;
    else
        YY1(i) = YY1(i) + 1900;
    end
end

YY2 = string(compose('%02d',YY1));
YY = str2double(YY2);

DD1 = mod(time,10000)-mod(time,100);
DD2 = DD1/100;
DD3 = string(compose('%02d',DD2));
DD = str2double(DD3);

MM1 = mod(time,1000000) - mod(time,10000);
MM2 = MM1/10000;
MM3 = string(compose('%02d',MM2));
MM = str2double(MM3);

hms(hms==-9) = nan;

ss1 = mod(hms,100);
ss2 = string(compose('%02d',ss1));
ss = str2double(ss2);

mm1 = mod(hms,10000)-mod(hms,100);
mm2 = mm1/100;
mm3 = string(compose('%02d',mm2));
mm = str2double(mm3);

hh1 = mod(hms,1000000) - mod(hms,10000);
hh2 = hh1/10000;
hh3 = string(compose('%02d',hh2));
hh = str2double(hh3);


time2 = datetime(YY,MM,DD,hh,mm,ss);

newCast = [1];

for i = 2:length(pB)
    if pB(i) < pB(i-1)
        disp(i);
        newCast = [newCast i];
    end
end

pgrid = nan(320,20);
tgrid = NaT(320,20);
chlagrid = NaN(320,20);

% First Case.
pgrid(1,1:(newCast(2)-newCast(1))) = pB(newCast(1):newCast(2)-1);
tgrid(1,1:(newCast(2)-newCast(1))) = time2(newCast(1):newCast(2)-1);
chlagrid(1,1:(newCast(2)-newCast(1))) = chla(newCast(1):newCast(2)-1);

% Then Loop.
for i = 2:319
    pgrid(i,1:(newCast(i+1)-newCast(i))) = pB(newCast(i):newCast(i+1)-1);
    tgrid(i,1:(newCast(i+1)-newCast(i))) = time2(newCast(i):newCast(i+1)-1);
    chlagrid(i,1:(newCast(i+1)-newCast(i))) = chla(newCast(i):newCast(i+1)-1);
end

% Final Case.
pgrid(320,1:(length(pB)+1-newCast(320))) = pB(newCast(320):length(pB));
tgrid(320,1:(length(time2)+1-newCast(320))) = time2(newCast(320):length(time2));
chlagrid(320,1:(length(chla)+1-newCast(320))) = chla(newCast(320):length(chla));

tgridDatenum = datenum(tgrid);

figure;
contourf(tgridDatenum,pgrid,chlagrid,'LineColor','auto');
set(gca,"YDir","reverse");
datetickzoom('x','yyyy','keeplimits');
colormap(flipud(cbrewer2("GnBu")));
c = colorbar;
c.Label.String = 'chl-a [ng/l]';
c.FontSize = 13;
ylim([1 18]);
zlim([0 500]);
yticks(1:1:18);
yticklabels(5:10:175);
ylabel("P [dbar]","FontSize",13); xlabel("Time",FontSize=13);
ax = gca;
ax.FontSize = 15;

%% Chl-a mean profile.

chlaProfile = nan(1,20);
f5 = nan(1,20);
f95 = nan(1,20);

for i = 1:20
    chlaProfile(i) = mean(chlagrid(pgrid==i),"omitnan");
    f5(i) = prctile(chlagrid(pgrid==i),5);
    f95(i) = prctile(chlagrid(pgrid==i),95);
end

% Toggle show y-label and title (for paper we don't need either)
displayYLabelAndTitle = false;

% Plot the mean profile of fluorescence with the mean and confidence
% interval.
figure;
plot(chlagrid(:,1),pgrid(:,1),'.',Color=[0.8 0.8 0.8],DisplayName="raw data");
hold on
plot(chlagrid(:,2:20),pgrid(:,2:20),'.',Color=[0.8 0.8 0.8],HandleVisibility='off');
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


%% Chl-a. Parameters vs Depth (assuming a lognormal distribution).


XN = nan(20,500);
for i = 1:20
    loc = find(pB == i);
    XN(i,1:length(loc)) = chla(loc);
end

figure;
histogram(XN(15,:));
XN(XN<=0) = nan;

for i = 1:20
    phat(i,:) = mle(XN(i,:),distribution="Lognormal");
end

depths = 5:10:195;

figure;
subplot(1,2,1)
plot(exp(phat(:,1)),depths);
title("$\mu^*$",Interpreter="latex"); set(gca,"YDir","reverse"); grid on
subplot(1,2,2)
plot(exp(phat(:,2)),depths);
title("$\sigma^*$",Interpreter="latex"); set(gca,"YDir","reverse"); grid on
sgtitle("L0 chl-$a$ parameters (1988-2021)",Interpreter="latex");

%% A-D Test.
if principle == true    
    
    % A-D
    tmpT = "-ad";
    
    % chla
    tmpX ="";
    tmp = importdata('data/L0/hplcChla_88-21_200.txt');
    [ax,~,~,pB,chla] = L0_helper(tmp,50,'ad');
    if showL0title == true
        sgtitle("L0 chl-$a$"+tmpX,"Interpreter","latex");
    end
    exportgraphics(ax,"figures/L0/bot/" + lp + "chla" + tmpT + ".png");
    clearvars -except tmpT startYear newCtd night logAxes lp seasonal;
    
    % chlb
    tmp = importdata('data/L0/chlb_88-21_200.txt');
    [ax,~,~] = L0_helper(tmp,50,'ad');
    sgtitle("L0: Chl $b$","Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "chlb" + tmpT + ".png");
    clearvars -except tmpT startYear newCtd night logAxes lp seasonal;
    
    % pc
    tmp = importdata('data/L0/pc_89-22_200.txt');
    [ax,~,~] = L0_helper(tmp,50,'ad');
    sgtitle("L0: Particulate Carbon","Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "pc" + tmpT + ".png");
    clearvars -except tmpT startYear newCtd night logAxes lp seasonal;
    
end
%% A-D Test: Seasonal.
if seasonal == true

    pVals = [];

    % WINTER
    tmpT = "-ad-01";

    % chla
    tmp = importdata('data/L0/hplcChla_88-21_200.txt');
    [ax,~,~,~,~,ad] = L0_helper(tmp,50,"ad",logAxes,1);
    sgtitle("L0: Chl $a$ (1988-2021)"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/"+lp+"chla" + tmpT + ".png");
    clear tmp ax; pVals = [pVals; ad([1 3])];

    % pc
    tmp = importdata('data/L0/pc_89-22_200.txt');
    [ax,~,~,~,~,ad] = L0_helper(tmp,50,"ad",logAxes,1);
    sgtitle("L0: Particulate Carbon"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "pc" + tmpT + ".png");
    clear tmp ax; pVals = [pVals; ad([1 3])];

    % SPRING
    tmpT = "-ad-02";

    % chla
    tmp = importdata('data/L0/hplcChla_88-21_200.txt');
    [ax,~,~,pB,chla,ad] = L0_helper(tmp,50,"ad",logAxes,2);
    sgtitle("L0: Chl $a$ (1988-2021)"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/"+lp+"chla" + tmpT + ".png");
    clear tmp ax; pVals = [pVals; ad([1 3])];

    % pc
    tmp = importdata('data/L0/pc_89-22_200.txt');
    [ax,~,~,~,~,ad] = L0_helper(tmp,50,"ad",logAxes,2);
    sgtitle("L0: Particulate Carbon"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "pc" + tmpT + ".png");
    clear tmp ax; pVals = [pVals; ad([1 3])];

    % SUMMER
    tmpT = "-ad-03";
    
    % chla
    tmp = importdata('data/L0/hplcChla_88-21_200.txt');
    [ax,~,~,~,~,ad] = L0_helper(tmp,50,"ad",logAxes,3);
    sgtitle("L0: Chl $a$ (1988-2021)"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/"+lp+"chla" + tmpT + ".png");
    clear tmp ax; pVals = [pVals; ad([1 3])];

    % pc
    tmp = importdata('data/L0/pc_89-22_200.txt');
    [ax,~,~,~,~,ad] = L0_helper(tmp,50,"ad",logAxes,3);
    sgtitle("L0: Particulate Carbon"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "pc" + tmpT + ".png");
    clear tmp ax; pVals = [pVals; ad([1 3])];

    % AUTUMN
    tmpT = "-ad-04";
    
    % chla
    tmp = importdata('data/L0/hplcChla_88-21_200.txt');
    [ax,~,~,~,~,ad] = L0_helper(tmp,50,"ad",logAxes,4);
    sgtitle("L0: Chl $a$ (1988-2021)"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/"+lp+"chla" + tmpT + ".png");
    clear tmp ax; pVals = [pVals; ad([1 3])];

    % pc
    tmp = importdata('data/L0/pc_89-22_200.txt');
    [ax,~,~,~,~,ad] = L0_helper(tmp,50,"ad",logAxes,4);
    sgtitle("L0: Particulate Carbon"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "pc" + tmpT + ".png");
    clear tmp ax; pVals = [pVals; ad([1 3])];

    season = [1 2 3 4];
    ax = figure;
    semilogy(season,pVals([1 18 35 52],:));
    lgd = legend("5","25",Location="south");
    lgd.Title.String = "pressure (dbar)";
    hold on; yline(0.05,"-","$\alpha$ = 0.05",Interpreter="latex",HandleVisibility="off"); yline(0.005,'--',"$\alpha$ = 0.005",Interpreter="latex",HandleVisibility="off");
    yline(0.1,'--',"$\alpha$ = 0.1",Interpreter="latex",HandleVisibility="off"); hold off; title("chl $a$ (A-D)",Interpreter="latex");
    set(gca,"XTick",1:1:4,"XTickLabel",["winter","spring","summer","autumn"]);
    xlabel("season"); ylabel("$p$-value",Interpreter="latex");
    exportgraphics(ax,"figures/L0/bot/synthSsnl/chla_ad.png");
end

%% A-D test: Chl-a, 2001-2021.
% Starting with cruise no. 131. This is to compare with the newer
% fluorometer on the CTD that entered use on this cruise.
if newCtd == true
    
    % A-D
    tmpT = "-ad";
    tmp = importdata('data/L0/hplcChla_01-22_200.txt');
    [ax,~,~] = L0_helper(tmp,50,'ad');
    sgtitle("L0: Chl $a$ (2001-2021)","Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "chla_01-21" + tmpT + ".png");
    clearvars -except tmpT startYear newCtd night logAxes lp;

else
    disp("Not analysing CRN 131 only data...");
end

%% A-D test: Chl-a, 2001-2022, night-time only.
if night == true
        
    tmpT = "";
    tmp = importdata('data/L0/hplcChla_01-22_200.txt');
    
    hms = char(string(tmp.data(:,3)));
    mdy = char(string(tmp.data(:,2)));
    
    % test = "0" + t2(1,1:end-1);
    for i = 1:length(tmp.data)
        if hms(i,end) == " "
            hms(i,:) = "0" + hms(i,1:end-1);
        end
        if mdy(i,end) == " "
            mdy(i,:) = "0" + mdy(i,1:end-1);
        end
    end
    
    Y = double("20" + mdy(:,5:6));
    M = double("" + mdy(:,1:2));
    D = double("" + mdy(:,3:4));
    h = double("" + hms(:,1:2));
    m = double("" + hms(:,3:4));
    s = double("" + hms(:,5:6));
    
    T = datetime(Y,M,D,h,m,s);
    T2 = datetime(Y,M,D);
    T3 = datenum(T);
    
    Lat = 22.75;
    Lon = -158;
    [SunRiseSet,~,~,~,~,~] = suncycle(Lat,Lon,T2);
    for i = 1:length(SunRiseSet)
        tmp = SunRiseSet(i,:) - 10;
        for j = 1:2
            if tmp(j) < 0
                tmp(j) = tmp(j) + 24;
            end
        end
        rs(i,:) = tmp;
    end
    rs2 = hours(rs);
    rs2.Format = 'hh:mm';
    
    sunriseTime = hours(rs2(:,1)');
    sunsetTime = hours(rs2(:,2)');
    
    dayFrac = rem(T3,1);
    castTime = dayFrac*24;
    
    castAtNight = nan(length(T),1);
    
    for i = 1:length(T)
        if castTime(i) < sunriseTime(i) || castTime(i) > sunsetTime(i)
            castAtNight(i) = 1;
        end
    end
    
    nightCastIDs = [];
    
    for i=1:length(T)
        if(~isnan(castAtNight(i)))
            disp(i);
            nightCastIDs = [nightCastIDs i];
        end
    end
    
    tmp = importdata('data/L0/hplcChla_01-22_200.txt').data(nightCastIDs,:);
    [ax,~,~] = L0_helper(tmp,50,"ks",logAxes);
    sgtitle("L0: Chl $a$ (2001-2021, NIGHT)","Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "chla_01-21_night" + tmpT + ".png");
    clearvars -except tmpT startYear newCtd night logAxes lp;

else
    disp("Not analysing night-time only data...");
end

%% A-D test: Chl-a. Start-Year Analysis.
% Here we move the start date of the analysis forward in time to see if
% the distribution of data has some dependence on time. It will only be
% checked if "startYear" is true.

if startYear==true
    
    % save p-values per year
    yearList = 1988:1:2016;
    pVals = [];
    
    tmpT = "-ad";
    
    % 2016-2022
    tmp = importdata('data/L0/hplcChla_16-22_200.txt');
    [ax,~,~,~,~,ad] = L0_helper(tmp,50,'ad');
    sgtitle("L0: Chl $a$ (2016-2022)"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "chla_16-22" + tmpT + ".png");
    clear ax tmp; pVals = [pVals; ad([1 3])];
    
    % 2015-2022
    tmp = importdata('data/L0/hplcChla_15-22_200.txt');
    [ax,~,~,~,~,ad] = L0_helper(tmp,50,'ad');
    sgtitle("L0: Chl $a$ (2015-2022)"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "chla_15-22" + tmpT + ".png");
    clear ax tmp; pVals = [pVals; ad([1 3])];
    
    % 2014-2022
    tmp = importdata('data/L0/hplcChla_14-22_200.txt');
    [ax,~,~,~,~,ad] = L0_helper(tmp,50,'ad');
    sgtitle("L0: Chl $a$ (2014-2022)"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "chla_14-22" + tmpT + ".png");
    clear ax tmp; pVals = [pVals; ad([1 3])];
    
    % 2013-2022
    tmp = importdata('data/L0/hplcChla_13-22_200.txt');
    [ax,~,~,~,~,ad] = L0_helper(tmp,50,'ad');
    sgtitle("L0: Chl $a$ (2013-2022)"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "chla_13-22" + tmpT + ".png");
    clear ax tmp; pVals = [pVals; ad([1 3])];
    
    % 2012-2022
    tmp = importdata('data/L0/hplcChla_12-22_200.txt');
    [ax,~,~,~,~,ad] = L0_helper(tmp,50,'ad');
    sgtitle("L0: Chl $a$ (2012-2022)"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "chla_12-22" + tmpT + ".png");
    clear ax tmp; pVals = [pVals; ad([1 3])];
    
    % 2011-2022
    tmp = importdata('data/L0/hplcChla_11-22_200.txt');
    [ax,~,~,~,~,ad] = L0_helper(tmp,50,'ad');
    sgtitle("L0: Chl $a$ (2011-2022)"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "chla_11-22" + tmpT + ".png");
    clear ax tmp; pVals = [pVals; ad([1 3])];
    
    % 2010-2022
    tmp = importdata('data/L0/hplcChla_10-22_200.txt');
    [ax,~,~,~,~,ad] = L0_helper(tmp,50,'ad');
    sgtitle("L0: Chl $a$ (2010-2022)"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "chla_10-22" + tmpT + ".png");
    clear ax tmp; pVals = [pVals; ad([1 3])];
    
    % 2009-2022
    tmp = importdata('data/L0/hplcChla_09-22_200.txt');
    [ax,~,~,~,~,ad] = L0_helper(tmp,50,'ad');
    sgtitle("L0: Chl $a$ (2009-2022)"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "chla_09-22" + tmpT + ".png");
    clear ax tmp; pVals = [pVals; ad([1 3])];
    
    % 2008-2022
    tmp = importdata('data/L0/hplcChla_08-22_200.txt');
    [ax,~,~,~,~,ad] = L0_helper(tmp,50,'ad');
    sgtitle("L0: Chl $a$ (2008-2022)"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "chla_08-22" + tmpT + ".png");
    clear ax tmp; pVals = [pVals; ad([1 3])];
    
    % 2007-2022
    tmp = importdata('data/L0/hplcChla_07-22_200.txt');
    [ax,~,~,~,~,ad] = L0_helper(tmp,50,'ad');
    sgtitle("L0: Chl $a$ (2007-2022)"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "chla_07-22" + tmpT + ".png");
    clear ax tmp; pVals = [pVals; ad([1 3])];
    
    % 2006-2022
    tmp = importdata('data/L0/hplcChla_06-22_200.txt');
    [ax,~,~,~,~,ad] = L0_helper(tmp,50,'ad');
    sgtitle("L0: Chl $a$ (2006-2022)"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "chla_06-22" + tmpT + ".png");
    clear ax tmp; pVals = [pVals; ad([1 3])];
    
    % 2005-2022
    tmp = importdata('data/L0/hplcChla_05-22_200.txt');
    [ax,~,~,~,~,ad] = L0_helper(tmp,50,'ad');
    sgtitle("L0: Chl $a$ (2005-2022)"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "chla_05-22" + tmpT + ".png");
    clear ax tmp; pVals = [pVals; ad([1 3])];
    
    % 2004-2022
    tmp = importdata('data/L0/hplcChla_04-22_200.txt');
    [ax,~,~,~,~,ad] = L0_helper(tmp,50,'ad');
    sgtitle("L0: Chl $a$ (2004-2022)"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "chla_04-22" + tmpT + ".png");
    clear ax tmp; pVals = [pVals; ad([1 3])];
    
    % 2003-2022
    tmp = importdata('data/L0/hplcChla_03-22_200.txt');
    [ax,~,~,~,~,ad] = L0_helper(tmp,50,'ad');
    sgtitle("L0: Chl $a$ (2003-2022)"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "chla_03-22" + tmpT + ".png");
    clear ax tmp; pVals = [pVals; ad([1 3])];
    
    % 2002-2022
    tmp = importdata('data/L0/hplcChla_02-22_200.txt');
    [ax,~,~,~,~,ad] = L0_helper(tmp,50,'ad');
    sgtitle("L0: Chl $a$ (2002-2022)"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "chla_02-22" + tmpT + ".png");
    clear ax tmp; pVals = [pVals; ad([1 3])];
    
    % 2001-2022 done above.
    tmp = importdata('data/L0/hplcChla_01-22_200.txt');
    [ax,~,~,~,~,ad] = L0_helper(tmp,50,'ad');
    sgtitle("L0: Chl $a$ (2001-2021)"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "chla_01-21" + tmpT + ".png");
    clear ax tmp; pVals = [pVals; ad([1 3])];

    % 2000-2022
    tmp = importdata('data/L0/hplcChla_00-22_200.txt');
    [ax,~,~,~,~,ad] = L0_helper(tmp,50,'ad');
    sgtitle("L0: Chl $a$ (2000-2022)"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "chla_00-22" + tmpT + ".png");
    clear ax tmp; pVals = [pVals; ad([1 3])];
    
    % 1999-2022
    tmp = importdata('data/L0/hplcChla_99-22_200.txt');
    [ax,~,~,~,~,ad] = L0_helper(tmp,50,'ad');
    sgtitle("L0: Chl $a$ (1999-2022)"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "chla_99-22" + tmpT + ".png");
    clear ax tmp; pVals = [pVals; ad([1 3])];
    
    % 1998-2022
    tmp = importdata('data/L0/hplcChla_98-22_200.txt');
    [ax,~,~,~,~,ad] = L0_helper(tmp,50,'ad');
    sgtitle("L0: Chl $a$ (1998-2022)"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "chla_98-22" + tmpT + ".png");
    clear ax tmp; pVals = [pVals; ad([1 3])];
    
    % 1997-2022
    tmp = importdata('data/L0/hplcChla_97-22_200.txt');
    [ax,~,~,~,~,ad] = L0_helper(tmp,50,'ad');
    sgtitle("L0: Chl $a$ (1997-2022)"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "chla_97-22" + tmpT + ".png");
    clear ax tmp; pVals = [pVals; ad([1 3])];
    
    % 1996-2022
    tmp = importdata('data/L0/hplcChla_96-22_200.txt');
    [ax,~,~,~,~,ad] = L0_helper(tmp,50,'ad');
    sgtitle("L0: Chl $a$ (1996-2022)"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "chla_96-22" + tmpT + ".png");
    clear ax tmp; pVals = [pVals; ad([1 3])];
    
    % 1995-2022
    tmp = importdata('data/L0/hplcChla_95-22_200.txt');
    [ax,~,~,~,~,ad] = L0_helper(tmp,50,'ad');
    sgtitle("L0: Chl $a$ (1995-2022)"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "chla_95-22" + tmpT + ".png");
    clear ax tmp; pVals = [pVals; ad([1 3])];
    
    % 1994-2022
    tmp = importdata('data/L0/hplcChla_94-22_200.txt');
    [ax,~,~,~,~,ad] = L0_helper(tmp,50,'ad');
    sgtitle("L0: Chl $a$ (1994-2022)"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "chla_94-22" + tmpT + ".png");
    clear ax tmp; pVals = [pVals; ad([1 3])];
    
    % 1993-2022
    tmp = importdata('data/L0/hplcChla_93-22_200.txt');
    [ax,~,~,~,~,ad] = L0_helper(tmp,50,'ad');
    sgtitle("L0: Chl $a$ (1993-2022)"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "chla_93-22" + tmpT + ".png");
    clear ax tmp; pVals = [pVals; ad([1 3])];
    
    % 1992-2022
    tmp = importdata('data/L0/hplcChla_92-22_200.txt');
    [ax,~,~,~,~,ad] = L0_helper(tmp,50,'ad');
    sgtitle("L0: Chl $a$ (1992-2022)"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "chla_92-22" + tmpT + ".png");
    clear ax tmp; pVals = [pVals; ad([1 3])];
    
    % 1991-2022
    tmp = importdata('data/L0/hplcChla_91-22_200.txt');
    [ax,~,~,~,~,ad] = L0_helper(tmp,50,'ad');
    sgtitle("L0: Chl $a$ (1991-2022)"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "chla_91-22" + tmpT + ".png");
    clear ax tmp; pVals = [pVals; ad([1 3])];
    
    % 1990-2022
    tmp = importdata('data/L0/hplcChla_90-22_200.txt');
    [ax,~,~,~,~,ad] = L0_helper(tmp,50,'ad');
    sgtitle("L0: Chl $a$ (1990-2022)"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "chla_90-22" + tmpT + ".png");
    clear ax tmp; pVals = [pVals; ad([1 3])];
    
    % 1989-2022
    tmp = importdata('data/L0/hplcChla_89-22_200.txt');
    [ax,~,~,~,~,ad] = L0_helper(tmp,50,'ad');
    sgtitle("L0: Chl $a$ (1989-2022)"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "chla_89-22" + tmpT + ".png");
    clear ax tmp; pVals = [pVals; ad([1 3])];

    % 1988-2022
    tmp = importdata('data/L0/hplcChla_88-22_200.txt');
    [ax,~,~,~,~,ad] = L0_helper(tmp,50,'ad');
    sgtitle("L0: Chl $a$ (1988-2022)"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "chla_89-22" + tmpT + ".png");
    clear ax tmp; pVals = [pVals; ad([1 3])];

    % start-year analysis synthesis plot
    ax = figure;
    subplot(1,2,1)
    semilogy(flip(yearList),pVals(1:29,:)); grid on;
    hold on; yline(0.05,"-","$\alpha$ = 0.05",Interpreter="latex"); yline(0.005,'--',"$\alpha$ = 0.005",Interpreter="latex");
    yline(0.1,'--',"$\alpha$ = 0.1",Interpreter="latex"); hold off; title("K-S"); xlim([1988 2016]);
    legend("5 dbar","25 dbar",Location="northwest"); xlabel("starting year"); ylabel("$p$-value",Interpreter="latex");
    subplot(1,2,2)
    semilogy(flip(yearList),pVals(30:58,:)); grid on;
    hold on; yline(0.05,"-","$\alpha$ = 0.05",Interpreter="latex"); yline(0.005,'--',"$\alpha$ = 0.005",Interpreter="latex");
    yline(0.1,'--',"$\alpha$ = 0.1",Interpreter="latex"); hold off; title("A-D"); xlim([1988 2016]);
    legend("5 dbar","25 dbar",Location="northwest"); xlabel("starting year"); ylabel("$p$-value",Interpreter="latex");
    sgtitle("does lognormality of surface chl $a$ depend on the start year of the analysis?","Interpreter","latex")
    exportgraphics(ax,"figures/L0/bot/startYear/surfaceChla.png");


else
    disp("Not analysing year by year...");
end