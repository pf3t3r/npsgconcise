% Calculate and save the depth of the DCM [1] across all HOT cruises.

clc; close all; clear;
% dcm = load("output/dcm.mat").dcm;
ctdData = load("datafiles\ctd_iso_ALL.mat").ctd;

%% Calculate the DCM depth
maxPcm = nan(329,1); minPcm = nan(329,1);
meanPcm = nan(329,1); stdPcm = nan(329,1);

for i = 1:329
    %fcm(i,:) = test.fcm;
    tmp = ctdData(i).pcm;
    if ~isempty(tmp)
        maxPcm(i) = max(tmp,[],'omitnan');
        minPcm(i) = min(tmp,[],'omitnan');
        meanPcm(i) = 2*round(mean(tmp/2,'omitnan'));
        stdPcm(i) = std(tmp,[],"omitnan");
    end
end

%% Save the DCM depth
save output/dcm.mat meanPcm minPcm maxPcm stdPcm -append;
clear i tmp;

% [1] the max, min, mean, and std of the DCM depth