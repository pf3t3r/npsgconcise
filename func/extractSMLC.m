function [idOut,pOut,XOut] = extractSMLC(id,p,X,pMaxMld)
%extractSMLC: extract the sub-mixed layer concentration given bottle ID id,
%pressure p, bottle concentration X, and pressure of the mixed layer depth
%pMaxMld.
% INPUTS
% id: bottle identifier
% p: pressure of that bottle (dbar)
% X: concentration of that bottle (various units)
% pMaxMld: pressure of the mixed layer depth (dbar)
% OUTPUTS
% idOut: bottle ID for a concentration within L2
% pOut: pressure (dbar) of a concentration within L2
% XOut: concentration (various units) of a quantity within L2

% Exclude data with no measurements (-9 here => NaN)
X(X==-9) = nan;
id = id(~isnan(X),:);
p = p(~isnan(X));
X1 = X(~isnan(X));

% Cruise number 'crn'
crn = str2num(id(:,1:3));

% Stop evaluating after crn = 329
for i = 1:length(crn)
    if crn(i) == 330
        stop = i;
        break
    elseif crn(i) > 330
        stop = i;
        break
    else
        stop = length(p) + 1;
    end
end

% Extract measurements below pMaxMld
L = stop-1; % No. of casts with chl measurements in cruise 1 - 329
tmpP_subML = nan(L,1);
tmpCRN_subML = nan(L,1);
tmpX_subML = nan(L,1);
botID = [];

% To find the measurements taken beneath pMaxMld, we populate the arrays
% just defined for cases where p > pMld...
for i = 1:L
    tmpMld = pMaxMld(crn(i));
    if p(i) > tmpMld
        tmpP_subML(i) = p(i);
        tmpCRN_subML(i) = crn(i);
        tmpX_subML(i) = X1(i);
        tmpStr = id(i,:);
        botID = [botID;tmpStr];
    end
end

% ...and remove nan values (which represent measurements above pMaxMld)
pOut = tmpP_subML(~isnan(tmpP_subML));
XOut = tmpX_subML(~isnan(tmpX_subML));
idOut = botID;

end