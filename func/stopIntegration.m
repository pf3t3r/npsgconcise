function [value, isterminal, direction] = stopIntegration(t2,C2)
    % this is an 'event' for use with the Andersson2 function.

    a = 1/1000;         % where we want to stop the integration
    value = C2(1) - a; 
    isterminal = 1;     % stop integraton
    direction = 0;      % all directions
end
