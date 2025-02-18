% Equation 2 from Andersson (2021). This equation describes the growth of
% phytoplankton. It assumes that growth comprises and deterministic
% constant component 'muk' and a stochastic component 'sigk*eta', where
% sigk is the magnitude of the stochastic variation and eta is the
% time-varying rate. For our analysis we introduce an additional mixing
% term 'M'. This can be constant or also time-varying. This equation is
% designed to be used with the ode45() solver.

function [dxdt,listeta] = andersson2(x,muk,sigk,M,timeVary)

% listeta was only used to verify that the random generation was working
% properly and is not necessary for the functioning of the model.

listeta = [];
if nargin < 5
    timeVary = false;
end
if nargin<4
    eta = normrnd(0,1);
    %eta = randi(1:2)-1.5;
    dxdt = (muk + sigk*eta)*x;
    listeta = [listeta eta];
else
    %eta = normrnd(0,1);
    eta = randi(1:2)-1.5;
    if timeVary == true
        dxdt = (muk + sigk*eta)*x - M*lognrnd(0,1);
    else
        dxdt = (muk + sigk*eta)*x - M;
    end
    listeta = [listeta eta];
end



end

