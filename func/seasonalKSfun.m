function [outputArg1,outputArg2] = seasonalKSfun(time,idStartMonth)
%UNTITLED work in progress

timeByMonth = discretize(month(time),12);

winter = [];
spring = [];
summer = [];
autumn = [];

for i=idStartMonth:329
    if timeByMonth(i) == 12 || timeByMonth(i) <= 2
        winter = [winter i];
    elseif timeByMonth(i) >= 3 && timeByMonth(i) <= 5
        spring = [spring i];
    elseif timeByMonth(i) >= 6 && timeByMonth(i) <= 8
        summer = [summer i];
    else
        autumn = [autumn i];
    end
end

% Try Eulerian first.
ksE_winter = zeros(5,depthMeasurements);
ksE_spring = zeros(5,depthMeasurements);
ksE_summer = zeros(5,depthMeasurements);
ksE_autumn = zeros(5,depthMeasurements);

% Try Lagrangian later.
ksL_winter = zeros(5,depthMeasurements);
ksL_spring = zeros(5,depthMeasurements);
ksL_summer = zeros(5,depthMeasurements);
ksL_autumn = zeros(5,depthMeasurements);

% Eulerian
for i = 1:depthMeasurements
    disp(i);
    tmp = chloro(i,winter);
    tmp(isnan(tmp) | tmp<=0) = [];
    [~,ksE_winter(:,i),~] = statsplot2(tmp,'noplot');
    tmp2 = chloro(i,spring);
    tmp2(isnan(tmp2) | tmp2<=0) = [];
    [~,ksE_spring(:,i),~] = statsplot2(tmp2,'noplot');
    tmp3 = chloro(i,summer);
    tmp3(isnan(tmp3) | tmp3<=0) = [];
    [~,ksE_summer(:,i),~] = statsplot2(tmp3,'noplot');
    tmp4 = chloro(i,autumn);
    tmp4(isnan(tmp4) | tmp4<=0) = [];
    [~,ksE_autumn(:,i),~] = statsplot2(tmp4,'noplot');
end
clear tmp;


% Lagrangian
for i = 1:depthMeasurements
    disp(i);
    tmp = chloroL(i,winter);
    tmp(isnan(tmp) | tmp<=0) = [];
    [~,ksL_winter(:,i),~] = statsplot2(tmp,'noplot');
    tmp2 = chloroL(i,spring);
    tmp2(isnan(tmp2) | tmp2<=0) = [];
    [~,ksL_spring(:,i),~] = statsplot2(tmp2,'noplot');
    tmp3 = chloroL(i,summer);
    tmp3(isnan(tmp3) | tmp3<=0) = [];
    [~,ksL_summer(:,i),~] = statsplot2(tmp3,'noplot');
    tmp4 = chloroL(i,autumn);
    tmp4(isnan(tmp4) | tmp4<=0) = [];
    [~,ksL_autumn(:,i),~] = statsplot2(tmp4,'noplot');
end


end