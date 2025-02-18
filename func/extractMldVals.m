function [idOut,pOut,XOut,crnOut] = extractMldVals(id,p,X,maxMld)
% extractMldVals() will extract vertical profiles of a bottle concentration
% within the mixed layer for each cruise at Station ALOHA. It will use an
% inputs the bottle ID, pressure, concentration of interest (e.g. Chl a), 
% and the maximum mixed layer depth 'maxMld' (defined per cruise). It will
% output vertical profiles within this mixed layer depth only. This output
% will be used for statistical tests on mixed layer dynamics.
% INPUTS
% id: bottle ID
% p: pressure [dbar]
% X: bottle concentration being measured [ng/l, ug/l]
% maxMld: maximum MLD of a cruise
% OUTPUTS
% idOut: bottle ID
% pOut: pressure [dbar]
% XOut: bottle concentration [ng/l, ug/l]

tmp = num2str(id);
bottleCRN = str2num(tmp(:,1:3));
clear tmp;

L = length(p);

% OPTION 2
tmpP = nan(1,L);
tmpCrn = nan(1,L);
tmpX = nan(1,L);
tmpId = nan(1,L);

for i = 1:L
    if bottleCRN(i) == 330
        stop = i;
        break
    elseif bottleCRN(i) > 330
        stop = i;
        break
    else
        stop = L+1;
    end
end

% disp(stop);

for i = 1:stop-1
    % OPTION 2: MAX MLD per cruise. This is what we will use.
    tmp = maxMld(bottleCRN(i));
    if p(i) < tmp
        tmpP(i) = p(i);
        tmpCrn(i) = bottleCRN(i);
        tmpX(i) = X(i);
        tmpId(i) = id(i);
    end
end

% OPTION 2
pOut = tmpP(~isnan(tmpP));
crnOut = tmpCrn(~isnan(tmpCrn));
XOut = tmpX(~isnan(tmpX));
idOut = tmpId(~isnan(tmpId));

end