clear; clc; addpath("func\"); close all; 

% Initial Conditions
x0 = 1;             % initial concentration
t0 = 0;             % start time
tF = 1000;          % finish time
h = 1;              % step size
time = t0:h:tF;     % time vector
r = 20;             % how many times we run each simulation

% Cases 1-2 run by default. Select which others are run by choosing
% true/false.
runCase3 = false;
runCase4 = false;
checkEtas = false;
runAnalytical = false;

%% Case 1: "Campbell"
% d[X]/dt = -k[X]
% In order to be consistent with later work, we set up 'k' as the 
% stochastic fluctuation in Andersson (2021). In other words, k =
% sigma*eta, where sigma is the magnitude of the stochastic fluctuation,
% and eta is a random number (normrnd(0,1)). Therefore the final equation
% used here is:
% d[X]/dt = + sigma*eta*[X]

sigma_k = 7.2;          % magnitude of the stochastic fluctuation

opt = odeset(Events=@stopIntegration);  % stop integration at 0.001*x0

% run the simulation of Case 1 'r' times
C1tmp = nan(length(time),r); t1 = nan(length(time),r);
for i = 1:r
    [tmpT,tmp] = ode45(@(t,x) andersson2(x,0,sigma_k),time,x0,opt);
    C1tmp(1:length(tmp),i) = tmp;
    t1(1:length(tmpT),i) = tmpT;
end

C1 = mean(C1tmp,2,"omitnan");
[tH1,tP1] = adtest(C1,"Distribution","logn");

%% Case 2: "Andersson". Deterministic + stochastic component to growth.
% Andersson (2021) proposes splitting the net growth rate k into a
% deterministic, constant mean component mu_k, and a stochastic fluctuation
% described by sigma_k*eta(t). Sigma_k is the magnitude of the stochastic
% fluctuation and eta(t) describes the time-dependency of the random
% fluctuations. In the original paper, a minus sign was used, but we here
% now use a plus sign.
% d[X]/dt = +(mu_k + sigma_k*eta(t))[X]

mu_k = 0.1;             % constant, deterministic component of growth

C2tmp = nan(length(time),r); t2 = nan(length(time),r);
for i = 1:r
    [tmpT,tmp] = ode45(@(t,x) andersson2(x,mu_k,sigma_k),time,x0,opt);
    C2tmp(1:length(tmp),i) = tmp;
    t2(1:length(tmpT),i) = tmpT;
end

if checkEtas == true
    for i = 1:r
        for k = 1:numel(t2(:,1))
            [~,etas(k,i)] = andersson2(C2tmp(k,i),mu_k,sigma_k);
        end
    end
end

C2 = mean(C2tmp,2,"omitnan");
[tH2,tP2] = adtest(C2,"Distribution","logn");

p = polyfit(time(C2>0),C2(C2>0),1);
f = polyval(p,time(C2>0));

%% Case 3. Add in the mixing as a constant.
% d[X]/dt = -(mu_k + sigma_k*eta(t))[X] + M.

if runCase3 == true
    C3tmp = nan(length(time),r);

    xTarget = 0.05;
    %mix = xTarget*mu_k;           % mixing term
    mix = 0.1;

    for i = 1:r
        [~,tmp] = ode45(@(t,x) andersson2(x,mu_k,sigma_k,mix),time,x0,opt);
        C3tmp(1:length(tmp),i) = tmp;
    end

    C3 = mean(C3tmp,2,"omitnan");

    [tH3,tP3] = adtest(C3,"Distribution","logn");
end

%% Case 4. Mixing varies with time.
% d[X]/dt = -(mu_k + sigma_k*eta(t))[X] + M(t).


if runCase4 == true
    C4tmp = nan(length(time),r);

    for i = 1:r
        [~,tmp] = ode45(@(t,x) andersson2(x,mu_k,sigma_k,mix,true),time,x0,opt);
        C4tmp(1:length(tmp),i) = tmp;
    end

    C4 = mean(C4tmp,2,"omitnan");
    [tH4,tP4] = adtest(C4,"Distribution","logn");
end

%% Compare Cases 1 - 4 visually.

% Plot also theoretical skewness/kurtosis for a lognormal distribution
sigTh = linspace(0,1,1000);
for i = 1:length(sigTh)
    skLogn(i) = (exp(sigTh(i)^2) + 2)*(sqrt(exp(sigTh(i)^2) - 1));
    kuLogn(i) = exp(4*sigTh(i)^2) + 2*exp(3*sigTh(i)^2) + 3*exp(2*sigTh(i)^2) - 3;
end

figure
subplot(3,3,[1 2 4 5 7 8])
semilogy(time,C1,LineStyle="-.",LineWidth=1,Color="#d95f02"); hold on
semilogy(time,C2,Color='#7570b3');
semilogy(time(f>0),f(f>0),Color="#000000",HandleVisibility="off");
% str = {'A simple plot','from 1 to 10'};
text(40,5,"numerical \mu = "+p(1));
if runCase3 == false && runCase4 == false
    title("C1: Campbell, C2: Andersson",FontSize=8);
    legend("C1","C2",Location="northwest");
elseif runCase3 == true
    semilogy(time,C3,Color='#e7298a');
    title("C1: Campbell, C2: Andersson, C3: A'son + mixing",FontSize=8);
end
if runCase4 == true
    semilogy(time,C4,Color='#66a61e');
    title("C1: Campbell, C2: Andersson, C3: A'son + mixing, C4: A'son + time-varying mixing",FontSize=8);
end
hold off; grid on
xlabel("time"); ylabel("concentration");
legend("C1","C2","C3","C4",Location="northwest");

subplot(3,3,3)
histogram(log(C1),DisplayName="C1",FaceColor="#d95f02"); hold on
h = histogram(log(C2),DisplayName="C2",FaceColor="#7570b3");
if runCase3 == true
histogram(log(C3),50,DisplayName="C3",FaceColor="#e7298a");
end
if runCase4 == true
histogram(log(C4),50,DisplayName="C4",FaceColor="#66a61e");
end
hold off
legend("C1","C2","C3","C4",Location="northeast");

subplot(3,3,[6 9])
plot(skLogn,kuLogn,DisplayName="Logn",Color="#1b9e77"); hold on
scatter(skewness(C1),kurtosis(C1),[],'MarkerEdgeColor','#000000','MarkerFaceColor','#d95f02');
scatter(skewness(C2),kurtosis(C2),[],'MarkerEdgeColor','#000000','MarkerFaceColor','#7570b3');
if runCase3 == true
    scatter(skewness(C3),kurtosis(C3),[],'MarkerEdgeColor','#000000','MarkerFaceColor','#e7298a');
end
if runCase4 == true
    scatter(skewness(C4),kurtosis(C4),[],'MarkerEdgeColor','#000000','MarkerFaceColor','#66a61e');
end
hold off
ylim([2 12]); xlim([0 3.5]);
legend('Log.','C1','C2','C3','C4',Location='southeast');

if runCase3 == true || runCase4 == true
    sgtitle("\mu_k = "+mu_k+", \sigma_k = "+sigma_k+", M = "+mix+"");
else
    sgtitle("\mu_k = "+mu_k+", \sigma_k = "+sigma_k+"");
end

%% Test basic solution of Case 1

if runAnalytical == true
    mu = 0; t = 0:(1000/100000):1000;
    sigma = 5e-4; eta = normrnd(0,1,length(t),1)';
    X1 = exp(-mu*t);             % case one
    X2 = exp(-(mu+sigma*eta).*t);  % case two
    
    figure
    plot(t,X1); hold on
    plot(t,X2); hold off
    disp(mean(X2));
    legend("case one","case two");
    title("\mu = " + mu + ", \sigma = " + sigma);
end

%% realistic mu study

muSig = [68.3 68.5 69 70 71 72 73 74 75 76 77 78 79 80 81];
numMu = [16.3508 10.4064 2.5199 4.3821 10.4134 0.28142 0.13845 2.3196 0.88711 0.16095 -0.0027239 0.17528 -0.002159 -0.0040504 -0.00095759];

figure
plot(muSig,numMu,'*-');
xlabel("\mu / \sigma");
ylabel("numerical \mu");
title("\mu / \sigma vs. numerical \mu");
subtitle("\sigma = 0.1");


%% theoretical lognormal
% x = 1;
% 
% binEdgesOfConc = h.BinEdges;
% 
% X = binEdgesOfConc;
% 
% P = (1./(X.*sqrt(2*pi*sigma_k.^2)) ) .* exp(-(log(X)-mu_k).^2 / 2*sigma_k.^2);

% figure;plot(X,P);