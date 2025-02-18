% Check how sensitive the Kolmogorov-Smirnov (K-S) test is to 
% - the number of significant digits in a dataset, and
% - the size of that dataset.

clc;close all;clear;

% Sample
sampleSize = [100,250,500,1000];

% Random Array for Each Sample Size
% ranNorm1000 = nan(1000,5);
% ranNorm1000(:,1) = randn(1000,1);

ranLogn1500 = nan(1500,5);
ranLogn1000 = nan(1000,5);
ranLogn500 = nan(500,5);
ranLogn250 = nan(250,5);
ranLogn100 = nan(100,5);

mu = 3.5; sd = 0.3;

ranLogn1500(:,1) = lognrnd(mu,sd,1500,1);
ranLogn1000(:,1) = lognrnd(mu,sd,1000,1);
ranLogn500(:,1) = lognrnd(mu,sd,500,1);
ranLogn250(:,1) = lognrnd(mu,sd,250,1);
ranLogn100(:,1) = lognrnd(mu,sd,100,1);

%%

figure;
subplot(2,1,1)
histogram(ranLogn1000);
subplot(2,1,2)
histogram(log(ranLogn1000));
%% Round

% 1500
ranLogn1500(:,2) = round(ranLogn1500(:,1),4,'significant');
ranLogn1500(:,3) = round(ranLogn1500(:,1),3,"significant");
ranLogn1500(:,4) = round(ranLogn1500(:,1),2,"significant");
ranLogn1500(:,5) = round(ranLogn1500(:,1),1,"significant");

% 1000
ranLogn1000(:,2) = round(ranLogn1000(:,1),4,"significant");
ranLogn1000(:,3) = round(ranLogn1000(:,1),3,"significant");
ranLogn1000(:,4) = round(ranLogn1000(:,1),2,"significant");
ranLogn1000(:,5) = round(ranLogn1000(:,1),1,"significant");

% 500
ranLogn500(:,2) = round(ranLogn500(:,1),4,"significant");
ranLogn500(:,3) = round(ranLogn500(:,1),3,"significant");
ranLogn500(:,4) = round(ranLogn500(:,1),2,"significant");
ranLogn500(:,5) = round(ranLogn500(:,1),1,"significant");

% 250
ranLogn250(:,2) = round(ranLogn250(:,1),4,"significant");
ranLogn250(:,3) = round(ranLogn250(:,1),3,"significant");
ranLogn250(:,4) = round(ranLogn250(:,1),2,"significant");
ranLogn250(:,5) = round(ranLogn250(:,1),1,"significant");

% 100
ranLogn100(:,2) = round(ranLogn100(:,1),4,"significant");
ranLogn100(:,3) = round(ranLogn100(:,1),3,"significant");
ranLogn100(:,4) = round(ranLogn100(:,1),2,"significant");
ranLogn100(:,5) = round(ranLogn100(:,1),1,"significant");


%% MLE

for i = 1:5
    mleK5(i,:) = mle(ranLogn1500(:,i),'distribution','Lognormal');
    mleK(i,:) = mle(ranLogn1000(:,i),'distribution','Lognormal');
    mle5(i,:) = mle(ranLogn500(:,i),'distribution','Lognormal');
    mle2(i,:) = mle(ranLogn250(:,i),'distribution','Lognormal');
    mle1(i,:) = mle(ranLogn100(:,i),'distribution','Lognormal');

    xCdfK5(i,:) = linspace(min(ranLogn1500(:,i)) - 2*std(ranLogn1500(:,i)), max(ranLogn1500(:,i)) + 2*std(ranLogn1500(:,i)),2000);
    xCdfK(i,:) = linspace(min(ranLogn1000(:,i)) - 2*std(ranLogn1000(:,i)), max(ranLogn1000(:,i)) + 2*std(ranLogn1000(:,i)),2000);
    xCdf5(i,:) = linspace(min(ranLogn500(:,i)) - 2*std(ranLogn500(:,i)), max(ranLogn500(:,i)) + 2*std(ranLogn500(:,i)),2000);
    xCdf2(i,:) = linspace(min(ranLogn250(:,i)) - 2*std(ranLogn250(:,i)), max(ranLogn250(:,i)) + 2*std(ranLogn250(:,i)),2000);
    xCdf1(i,:) = linspace(min(ranLogn100(:,i)) - 2*std(ranLogn100(:,i)), max(ranLogn100(:,i)) + 2*std(ranLogn100(:,i)),2000);

    yCdfLognK5(i,:) = cdf('logn', xCdfK5(i,:), mleK5(i,1), mleK5(i,2));
    yCdfLognK(i,:) = cdf('logn', xCdfK(i,:), mleK(i,1), mleK(i,2));
    yCdfLogn5(i,:) = cdf('logn', xCdf5(i,:), mle5(i,1), mle5(i,2));
    yCdfLogn2(i,:) = cdf('logn', xCdf2(i,:), mle2(i,1), mle2(i,2));
    yCdfLogn1(i,:) = cdf('logn', xCdf1(i,:), mle1(i,1), mle1(i,2));

end

%% KS

for i = 1:5
    [Hl_k5(i),KS_k5(i)] = kstest(ranLogn1500(:,i),[xCdfK5(i,:)' yCdfLognK5(i,:)']);
    [Hl_k(i),KS_k(i)] = kstest(ranLogn1000(:,i),[xCdfK(i,:)' yCdfLognK(i,:)']);
    [Hl_5(i),KS_5(i)] = kstest(ranLogn500(:,i),[xCdf5(i,:)' yCdfLogn5(i,:)']);
    [Hl_2(i),KS_2(i)] = kstest(ranLogn250(:,i),[xCdf2(i,:)' yCdfLogn2(i,:)']);
    [Hl_1(i),KS_1(i)] = kstest(ranLogn100(:,i),[xCdf1(i,:)' yCdfLogn1(i,:)']);
end

%%
noSig = [5 4 3 2 1];

ax1 = figure;
plot(noSig,KS_k5,'DisplayName','1500');
hold on
plot(noSig,KS_k,'DisplayName','1000');
plot(noSig,KS_5,'DisplayName','500');
plot(noSig,KS_2,'DisplayName','250');
plot(noSig,KS_1,'DisplayName','100');
hold off
grid on;
xlabel('No. of significant digits'); xlim([1 4]);
ylabel('p-value'); ylim([0 1]);
title('Precision Sensitivity for KS (random lognormal: mu = 3.5, sigma = 0.3)');
lgd = legend('Location','best');
lgd.Title.String = 'No. of samples';
exportgraphics(ax1,'figures/prec/testPrecS_3p5.png');

clear i;