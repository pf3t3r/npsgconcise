function [mleParams,pKs,pLil,pAd,pSw] = statsplotComp(x)
%
% Calculates Maximum Likelihood Estimates (MLE) of the parameters from 
% Normal, Lognormal, Gamma, Weibull and Exponential distributions; and 
% performs a Kolmogorov Smirnov Test (K-S) according to the output
% parameters. This code additionally performs the Lilliefors-corrected K-S
% Test (Lil), an Anderson-Darling Test (A-D), and a Shapiro-Wilks Test
% (S-W) on the data. The resultant p-values may be compared to each other.
%
% INPUT
%   x: the vector of unknown distribution.
% OUTPUT
%   mleParams: a 5*2 matrix of the parameters of the MLE. Rows are,
%       respectively, Normal, Lognormal, Gamma, Weibull and Exponential.
%       Last one has only one parameter assigned. Required only for the
%       Kolmogorov-Smirnov Test (K-S).
%   pKs: a 5*1 vector containing p-values from the Kolmogoroff-Smirnov test
%       for data versus MLE cdf.
%   pLil: a 5*1 vector containing p-values from the Lilliefors-corrected
%       K-S test on data.
%   pAd: a 5*1 vector containing p-values from the Anderson-Darling test
%       on data.
%   pSw: a 2*1 vector containing p-values from the Shapiro-Wilks test on 
%       data (normal and lognormal distribution only).
% Author: B. Barone, P. Farrell.

% Initialize Output values
mleParams = nan*ones(5,2);
pKs = nan*ones(5,1); pLil = nan*ones(5,1); pAd = nan*ones(5,1); pSw = nan*ones(2,1);
[mleParams(1,:),~] = mle(x,'distribution','norm','Alpha',0.32);
if sum(x<=0) == 0
    [mleParams(2,:),~] = mle(x,'distribution','logn','Alpha',0.32);
    [mleParams(3,:),~] = mle(x,'distribution','wbl','Alpha',0.32);
    [mleParams(4,:),~] = mle(x,'distribution','gamma','Alpha',0.32);
    [mleParams(5,1),~] = mle(x,'distribution','exp','Alpha',0.32);
else 
    mleParams(2,:) = [nan nan];
    disp('Can''t compute MLE parameters for Lognormal, Weibull, Gamma and Exponential distribution. Negative or zero values in input vector.')
end

% Compute cdfs for MLE distributions
x_cdf = linspace(min(x) - 2*std(x), max(x) + 2*std(x), 2000);
y_cdf_norm = cdf('norm', x_cdf, mleParams(1,1), mleParams(1,2));
y_cdf_logn = cdf('logn', x_cdf, mleParams(2,1), mleParams(2,2));
y_cdf_wbl = cdf('wbl', x_cdf, mleParams(3,1), mleParams(3,2));
y_cdf_gam = cdf('gamma', x_cdf, mleParams(4,1), mleParams(4,2));
y_cdf_exp = cdf('exp', x_cdf, mleParams(5,1));

% K-S test for each MLE distribution
[~,pKs(1)] = kstest(x,[x_cdf' y_cdf_norm']);
if sum(x<=0) == 0
    [~,pKs(2)] = kstest(x,[x_cdf' y_cdf_logn']);
    [~,pKs(3)] = kstest(x,[x_cdf' y_cdf_wbl']);
    [~,pKs(4)] = kstest(x,[x_cdf' y_cdf_gam']);
    [~,pKs(5)] = kstest(x,[x_cdf' y_cdf_exp']);
end

% NOTE that the MLE is not required for the following tests.

% Lil test on data
[~,pLil(1)] = lillietest(x,"MCTol",1e-2);
if sum(x<=0) == 0
    [~,pLil(2)] = lillietest(log(x),"MCTol",1e-2);
    [~,pLil(3)] = lillietest(log(x),"Distr","ev","MCTol",1e-2);
    [~,pLil(4)] = lillietest(x,MCTol=1e-2);
    [~,pLil(5)] = lillietest(x,MCTol=1e-2);
end

% create gamma PDO
pd2 = fitdist(x',"Gamma");

% A-D test on data
[~,pAd(1)] = adtest(x,MCTol=1e-2);
if sum(x<=0) == 0
    [~,pAd(2)] = adtest(x,"Distribution","logn",MCTol=1e-2);
    [~,pAd(3)] = adtest(x,"Distribution","weibull",MCTol=1e-2);
    [~,pAd(4)] = adtest(x,Distribution=pd2,MCTol=0.05);
    [~,pAd(5)] = adtest(x,MCTol=1e-2);
end

% S-W Test on data
[~,pSw(1)] = swtest(x);
[~,pSw(2)] = swtest(log(x));

% % Negative Log-Likelihood for each MLE distribution on data
% nll(1) = normlike(mleParams(1,:),x);
% nll(2) = lognlike(mleParams(2,:),x);
% nll(3) = wbllike(mleParams(3,:),x);
% nll(4) = gamlike(mleParams(4,:),x);
% nll(5) = explike(mleParams(5,1),x);

% % Plotting Results
% if strcmp(nc,'noplot') == 0
%     n_x = length(x);
%     x_pdf = linspace(min(x)-2*std(x),max(x)+2*std(x),2000);
%     y_pdf_norm = pdf('norm',x_pdf,mleParams(1,1),mleParams(1,2));
%     y_pdf_logn = pdf('logn',x_pdf,mleParams(2,1),mleParams(2,2));
%     y_pdf_wbl = pdf('wbl',x_pdf,mleParams(3,1),mleParams(3,2));
%     y_pdf_gam = pdf('gamma',x_pdf,mleParams(4,1),mleParams(4,2));
%     y_pdf_exp = pdf('exp',x_pdf,mleParams(5,1));
%     sp1 = subplot(2,3,1);
%     [n_h,x_h] = hist(x,nc);
%     hist(x,nc)
%     step = x_h(end)-x_h(end-1);
%     hold on
%     plot(x_pdf,y_pdf_norm*step*n_x,'r','linewidth',2)
%     plot(x_pdf,y_pdf_logn*step*n_x,'g','linewidth',2)
%     plot(x_pdf,y_pdf_wbl*step*n_x,'c','linewidth',2)
%     plot(x_pdf,y_pdf_gam*step*n_x,'m','linewidth',2)
%     hold off
%     title('\bf Histogram & MLE pdf')
%     legend('data','norm','logn','weib','gam')
%     ylabel('\bf n. of values')
%     xlabel('\bf value')
%     xlim([min(x)-std(x) max(x)+std(x)])
%     sp2 = subplot(2,3,2);
%     [ecdf_f,ecdf_x] = ecdf(x);
%     plot(ecdf_x,1-ecdf_f,'linewidth',2)
%     hold on
%     plot(x_cdf,1-y_cdf_norm,'r',x_cdf,1-y_cdf_logn,'g',x_cdf,1-y_cdf_wbl,'c',x_cdf,1-y_cdf_gam,'m')
%     hold off
%     xlim([min(x) max(x)])
%     title('\bf Empirical & MLE cdfs')
%     sp3 = subplot(2,3,3);
%     plot(ecdf_x,1-ecdf_f,'linewidth',2)
%     hold on
%     plot(x_cdf,1-y_cdf_norm,'r',x_cdf,1-y_cdf_logn,'g',x_cdf,1-y_cdf_wbl,'c',x_cdf,1-y_cdf_exp,'y')
%     hold off
%     set(gca,'xscale','log','yscale','log')
%     xlim([mean(x)+std(x) max(x)+std(x)])
%     ylim([0.00001 0.1])
%     legend('data','norm','logn','weib','exp','location','southwest')
%     title('\bf log-log Empirical & MLE cdfs')
%     sp4 = subplot(2,3,4);
%     hpp = probplot('normal',x);
%     set(hpp(1),'Markersize',4,'marker','o','color','k')
%     set(hpp(2),'Linestyle','-')
%     annotation('textbox',[0.1582    0.4433    0.1066    0.0417],'string',['KS-p: ' num2str(pKs(1),3)],'Edgecolor','none');
%     set(sp1,'Position',[0.1823    0.5838    0.2097    0.3412])
%     %set(sp2,'Position',[0.4108    0.5838    0.2134    0.3412])
%     set(sp3,'Position',[0.6502    0.5838    0.2097    0.3412])
%     set(sp4,'Position',[0.1580    0.1387    0.2097    0.3412])
%     if sum(x<=0) == 0
%         sp5 = subplot(2,3,5);
%         hpp = probplot('lognormal',x);
%         set(hpp(1),'Markersize',4,'marker','o','color','k')
%         set(hpp(2),'Linestyle','-')
%         ylabel('')
%         sp6 = subplot(2,3,6);
%         hpp = probplot('weibull',x);
%         set(hpp(1),'Markersize',4,'marker','o','color','k')
%         set(hpp(2),'Linestyle','-')
%         ylabel('')
%         set(sp5,'Position',[0.4145    0.1400    0.2097    0.3412])
%         set(sp6,'Position',[0.6693    0.1413    0.2097    0.3412])
%         annotation('textbox',[0.4168    0.4407    0.0845    0.0417],'string',['KS-p: ' num2str(pKs(2),3)],'Edgecolor','none');
%         annotation('textbox',[0.6651    0.4433    0.0953    0.0417],'string',['KS-p: ' num2str(pKs(3),3)],'Edgecolor','none');
%     end
end