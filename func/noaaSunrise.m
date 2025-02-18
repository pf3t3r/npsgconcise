function [sunrise,sunset] = noaaSunrise(dayOfYear,hour,year,lat,lon,timezone)
% Calculate the sunrise and sunset times for days input
% INPUT
% lon: longitude in degrees (positive to east of prime meridian, ALOHA =
% -158)
% timezone: hours from UTC (HT = -10)
% OUTPUT
% sunrise [UTC, min]
% sunset [UTC, min]

if mod(year,4) == 0
    yrLength = 366;
else
    yrLength = 365;
end

% fractional year 'gamma' [rad]
gamma = ((2*pi)/(yrLength)) * (dayOfYear - 1 + (hour - 12)/24);

% eqTime [min]
eqTime = 229.18*(0.000075 + 0.001868*cos(gamma) - 0.032077*sin(gamma) - ...
    0.014615*cos(2*gamma) - 0.040849*sin(2*gamma));
decl = 0.006918 - 0.399912*cos(gamma) + 0.070257*sin(gamma) - ...
    0.006758*cos(2*gamma) + 0.000907*sin(2*gamma) - 0.002697*cos(3*gamma) + ...
    0.00148*sin(3*gamma);
% timeOffset = eqTime + 4*lon - 60*timezone;
% tst = hour*60 + min + (sec/60) + timeOffset;
% solar hour angle [deg]
% ha = (tst/4) - 180;

% at sunrise/sunset, zenith angle = 90.833 deg
% convert all to rad
zen = deg2rad(90.833);
lat = deg2rad(lat);
decl = deg2rad(decl);

% solar hour angle for sunrise/sunset [deg]
% positive = sunrise, negative = sunset
ha = [acos( (cos(zen)/(cos(lat)*cos(decl))) - tan(lat)*tan(decl) ) ...
    -acos( (cos(zen)/(cos(lat)*cos(decl))) - tan(lat)*tan(decl) )];

sunrise = 720 - 4*(lon + rad2deg(ha(1))) - eqTime;
sunset = 720 - 4*(lon + rad2deg(ha(2))) - eqTime;

sunrise = minutes(sunrise);
sunset = minutes(sunset);
sunrise.Format = "hh:mm";
sunset.Format = 'hh:mm';

sunrise = sunrise + hours(timezone);
sunset = sunset + hours(timezone);

end