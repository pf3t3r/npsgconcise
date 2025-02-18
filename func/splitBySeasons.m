function [winter,spring,summer,autumn] = splitBySeasons(X,n)
%splitBySeasons split time array X of length n by seasons
% INPUT: time vector X of length n
% OUTPUT: save ID of vector X that corresponds to winter, spring, summer,
% or autumn respectively.

timeByMonth = month(X);

winter = nan; spring = nan; summer = nan; autumn = nan;

for i = 1:n
    if timeByMonth(i) <= 3
        winter = 1;
    elseif timeByMonth(i) >= 4 && timeByMonth(i) <= 6
        spring = 1;
    elseif timeByMonth(i) >= 7 && timeByMonth(i) <= 9
        summer = 1;
    else
        autumn = 1;
    end
end

end