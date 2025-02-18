function [idF,sigF,xF,crn] = densityConcentrations(id, p, X, T, Sp, pMaxMld)
% This function determines where a bottle concentration has coincident
% salinity and temperature measurements and finds the conservative
% temperature, absolute salinity, and density (referenced to surface
% pressure).

% INPUTS
% id: bottle ID
% p: pressure
% X: bottle concentration
% T: temperature (C)
% Sp: practical salinity (g/kg)
% OUTPUTS
% idF: bottle ID of filtered data
% sigF: in-situ density (kg/m^3), ref. pressure = 0 dbar,
% xF: bottle concentration of filtered data

X(X==-9) = nan;

id = id(~isnan(X),:);
p = p(~isnan(X));
X1 = X(~isnan(X));
T = T(~isnan(X));
Sp = Sp(~isnan(X));

% pmaxmld

crn = str2num(id(:,1:3));

for i = 1:length(crn)
    if crn(i) == 330
        stop = i;
        break
    else
        stop = length(p) + 1;
    end
end

%disp(stop);

% Extract meas below pMaxMld
L = stop-1; % No. of casts with chl measurements in cruise 1 - 329

tmpP_subML = nan(L,1);
tmpCRN_subML = nan(L,1);
tmpChl_subML = nan(L,1);
% tmpCrnStr_subML = nan(L,1);
tmpT = nan(L,1); tmpS = nan(L,1);
botID = [];

% To find the measurements taken beneath pMaxMld, we save them to a new
% file...
for i = 1:L
    %crn(i)
    % with MAX MLD per cruise
    tmpMld = pMaxMld(crn(i));
    if p(i) > tmpMld
        tmpP_subML(i) = p(i);
        tmpCRN_subML(i) = crn(i);
        tmpChl_subML(i) = X1(i);
        tmpStr = id(i,:);
        tmpT(i) = T(i);
        tmpS(i) = Sp(i);
        botID = [botID;tmpStr];
    end
end

% ...and remove nan values (which here represent measurements above
% pMaxMld)
pSubML = tmpP_subML(~isnan(tmpP_subML));
% crnSubML = tmpCRN_subML(~isnan(tmpCRN_subML));
chlSubML = tmpChl_subML(~isnan(tmpChl_subML));
tSubML = tmpT(~isnan(tmpT));
sSubML = tmpS(~isnan(tmpS));

% Calculate the sub-ML Density 'sig0_sM'
disp(length(sSubML));
disp(length(tSubML));
disp(length(pSubML));
SA_subML = gsw_SA_from_SP(sSubML,pSubML,158,22.75);
CT_subML = gsw_CT_from_t(SA_subML,tSubML,pSubML);

sig0_sM = gsw_sigma0(SA_subML,CT_subML);

% Bin the sub-ML density -> 'sig0_sMr'

sig0_sMr = round(sig0_sM,3,"significant");

% Display Histogram of the binned sub-ML density 'sig0_sMr' 

% Edges = (no. of unique values/bins + missing bins) - 1
% uniH = length(unique(sig0_sMr));
edgesH = length(min(sig0_sMr):0.1:max(sig0_sMr)) - 1;

figure;
histogram(sig0_sMr,edgesH);

% Remove densest waters -> continuous, binned sub-ML density 'sig0_sMrf'

% This is justified because (1) it will give a continuous 'depth' profile
% of density, and (2) 99.4% of the measurements are contained within the
% continuosu portion anyway.

sig0_sMrf = sig0_sMr(sig0_sMr<=25.8);

% Confirm that filtering has been done correctly
figure; histogram(sig0_sMrf,'BinWidth',0.1,'BinLimits',[min(sig0_sMrf) max(sig0_sMrf)]);

% Extract Chl-a at the filtered, binned densities found above

chlaSubMLf = chlSubML(sig0_sMr<=25.8);
botIDf = botID(sig0_sMr<=25.8,:);             % this is needed for ksOfLagrangian()
% crnSubMLf = crnSubML(sig0_sMr<=25.8);

idF = botIDf;
sigF = sig0_sMrf;
xF = chlaSubMLf;
end