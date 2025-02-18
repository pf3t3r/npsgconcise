function [pb5,pb10,X_out,n5,n10] = cleanAndBin(p,X,botid)
%cleanAndBin
% p = pressure [dbar]
% X = concentration of substance you want to measure
% botid = bottle ID

% Remove bottles that are too close to the surface (< 2.5 dbar)
idRm = p > 2.5;
p = p(idRm);
X = X(idRm);
botid = botid(idRm);

% Remove bottles where concentration of X = 0
idZero = X == 0;
p = p(~idZero);
X = X(~idZero);
botid = botid(~idZero);

% Save cruise number (CRN) of each bottle - needed below
tmp = num2str(botid);
bottleCRN = str2num(tmp(:,1:3));
clear tmp;

% Remove bottles from cruises 330 on (b/c fluorescence analysis not done)
for i = 1:length(p)
    if bottleCRN(i) > 329
        id329 = i - 1;
        break;
    else
        id329 = length(p);
    end
end

p = p(1:id329);
X = X(1:id329);
clear idRm idZero id329 i;

pb5 = discretize(p,0:5:200);
pb10 = discretize(p,0:10:200);
n5 = max(pb5);
n10 = max(pb10);

X_out = X;

end