close all; clc; clear;

dcm = load("output/dcm.mat").dcm; % pDcm + sigmaDcm (all casts, crn 1-329)

crn = dcm(:,1);
pcm = dcm(:,3);
dcmByCrn = nan(329,40);

k = 1;
while k < 329
    tmp = []; 
    for i = 1:10199
        if crn(i) == k
            tmp = [tmp pcm(i)];
        end
        
    end
    dcmByCrn(k,1:length(tmp)) = tmp;
    k = k + 1;
end

meanDcm = mean(dcmByCrn,2,"omitnan");

timeSeriesMeanDcm = mean(meanDcm,"omitnan");
timeSeriesPrctl = prctile(meanDcm,[5 95]);

save datafiles/meanDcm.mat meanDcm;
save datafiles/timeSeriesMeanDcm.mat timeSeriesMeanDcm;
save datafiles/timeSeriesPrctl.mat timeSeriesPrctl;