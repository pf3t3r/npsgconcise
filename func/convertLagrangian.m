function [lagrangianData] = convertLagrangian(X,o,n)
%convertLagrangian converts data in depth/pressure coordinates into
%Lagrangian coordinates, i.e. coordinates that are centered on the Deep
% Chlorophyll Maximum (DCM).
% INPUT: input vector 'X' in depth/pressure coordinates; amount that the 
% vector must be moved in the y direction in order to centre the DCM 'o';
% length of dataset 'n'.
% OUTPUT: lagrangianData

for k = 1:length(n)
    lagrangianData(:,k) = circshift(X,o);
    if o > -1 && o < 40
        lagrangianData(1:o,k) = NaN;
    elseif o == -1
        lagrangianData(end,k) = NaN;
    elseif o < -1 && o > -40
        lagrangianData((end+o):end,k) = NaN;
    elseif abs(o) > 40
        lagrangianData(:,k) = NaN;
    end
end