function [r1,r2,r3,r4] = quickRound(input)
%quickRound Rounds things quickly
r1 = round(input(:,1),4,'significant');
r2 = round(input(:,1),3,'significant');
r3 = round(input(:,1),2,'significant');
r4 = round(input(:,1),1,'significant'); 
end