function [MLEp,KSp,nll] = statsplot2(x,nc)
%
%
% function [MLEp,KSp,nll] = statsplot2(x)
%          [MLEp,KSp,nll] = statsplot2(x,nc)
%          [MLEp,KSp,nll] = statsplot2(x,'noplot')
% 
% Calculates Maximum Likelihood Estimates of the parameters from Normal,
% Lognormal, Gamma, Weibull and Exponential distributions; performs a 
% Kolmogorov Smirnov test on data for each MLE distribution; calculates 
% negative log-likelihood for each distribution to produce data; plots 
% histograms of the input vector, its empirical cdf and p-p plots for three
% distributions.
%
% INPUT
%   x: the vector of unknown distribution.
%   nc(OPTIONAL): number of classes for histogram plot (default = 20).
%       nc can also be assigned the value 'noplot' not to display results
%       graphically.           
% OUTPUT
%   MLEp: a 5*2 matrix of the parameters of the MLE. Rows are,
%       respectively, Normal, Lognormal, Gamma, Weibull and Exponential.
%       Last one has only one parameter assigned.
%   KSp: a 5*1 vector containing p-values from the Kolmogoroff-Smirnov test
%       for data versus MLE cdf.
%   nll: a 5*1 matrix with the negative log-likelihood for MLE distribution
%       to produce observed data.
%
% Author: B.Barone
%

% Assign nc if not user-defined
if nargin < 2
    nc = 20;
end
% Initialize Output values
MLEp = NaN*ones(5,2);
KSp = NaN*ones(5,1);
nll = NaN*ones(5,1);
% Finding MLE coefficients
[MLEp(1,:),~] = mle(x,'distribution','norm','Alpha',0.32);
if sum(x<=0) == 0
    [MLEp(2,:),~] = mle(x,'distribution','logn','Alpha',0.32);
    [MLEp(3,:),~] = mle(x,'distribution','wbl','Alpha',0.32);
    [MLEp(4,:),~] = mle(x,'distribution','gamma','Alpha',0.32);
    [MLEp(5,1),~] = mle(x,'distribution','exp','Alpha',0.32);
else MLEp(2,:) = [NaN NaN];
    display('Can''t compute MLE parameters for Lognormal, Weibull, Gamma and Exponential distribution. Negative or zero values in input vector.')
end

% Compute cdfs for MLE distributions
x_cdf = linspace(min(x)-2*std(x),max(x)+2*std(x),2000);
y_cdf_norm = cdf('norm',x_cdf,MLEp(1,1),MLEp(1,2));
y_cdf_logn = cdf('logn',x_cdf,MLEp(2,1),MLEp(2,2));
y_cdf_wbl = cdf('wbl',x_cdf,MLEp(3,1),MLEp(3,2));
y_cdf_gam = cdf('gamma',x_cdf,MLEp(4,1),MLEp(4,2));
y_cdf_exp = cdf('exp',x_cdf,MLEp(5,1));
% Kolmogorov Smirnov test for each MLE distribution
[Hx,KSp(1)] = kstest(x,[x_cdf' y_cdf_norm']);
if sum(x<=0) == 0
    [Hl,KSp(2)] = kstest(x,[x_cdf' y_cdf_logn']);
    [Hw,KSp(3)] = kstest(x,[x_cdf' y_cdf_wbl']);
    [Hg,KSp(4)] = kstest(x,[x_cdf' y_cdf_gam']);
    [He,KSp(5)] = kstest(x,[x_cdf' y_cdf_exp']);
end
% Negative Log-Likelihood for each MLE distribution on data
nll(1) = normlike(MLEp(1,:),x);
nll(2) = lognlike(MLEp(2,:),x);
nll(3) = wbllike(MLEp(3,:),x);
nll(4) = gamlike(MLEp(4,:),x);
nll(5) = explike(MLEp(5,1),x);
% Plotting Results
if strcmp(nc,'noplot') == 0
    n_x = length(x);
    x_pdf = linspace(min(x)-2*std(x),max(x)+2*std(x),2000);
    y_pdf_norm = pdf('norm',x_pdf,MLEp(1,1),MLEp(1,2));
    y_pdf_logn = pdf('logn',x_pdf,MLEp(2,1),MLEp(2,2));
    y_pdf_wbl = pdf('wbl',x_pdf,MLEp(3,1),MLEp(3,2));
    y_pdf_gam = pdf('gamma',x_pdf,MLEp(4,1),MLEp(4,2));
    y_pdf_exp = pdf('exp',x_pdf,MLEp(5,1));
    sp1 = subplot(2,3,1);
    [n_h,x_h] = hist(x,nc);
    hist(x,nc)
    step = x_h(end)-x_h(end-1);
    hold on
    plot(x_pdf,y_pdf_norm*step*n_x,'r','linewidth',2)
    plot(x_pdf,y_pdf_logn*step*n_x,'g','linewidth',2)
    plot(x_pdf,y_pdf_wbl*step*n_x,'c','linewidth',2)
    plot(x_pdf,y_pdf_gam*step*n_x,'m','linewidth',2)
    hold off
    title('\bf Histogram & MLE pdf')
    legend('data','norm','logn','weib','gam')
    ylabel('\bf n. of values')
    xlabel('\bf value')
    xlim([min(x)-std(x) max(x)+std(x)])
    sp2 = subplot(2,3,2);
    [ecdf_f,ecdf_x] = ecdf(x);
    plot(ecdf_x,1-ecdf_f,'linewidth',2)
    hold on
    plot(x_cdf,1-y_cdf_norm,'r',x_cdf,1-y_cdf_logn,'g',x_cdf,1-y_cdf_wbl,'c',x_cdf,1-y_cdf_gam,'m')
    hold off
    xlim([min(x) max(x)])
    title('\bf Empirical & MLE cdfs')
    sp3 = subplot(2,3,3);
    plot(ecdf_x,1-ecdf_f,'linewidth',2)
    hold on
    plot(x_cdf,1-y_cdf_norm,'r',x_cdf,1-y_cdf_logn,'g',x_cdf,1-y_cdf_wbl,'c',x_cdf,1-y_cdf_exp,'y')
    hold off
    set(gca,'xscale','log','yscale','log')
    xlim([mean(x)+std(x) max(x)+std(x)])
    ylim([0.00001 0.1])
    legend('data','norm','logn','weib','exp','location','southwest')
    title('\bf log-log Empirical & MLE cdfs')
    sp4 = subplot(2,3,4);
    hpp = probplot('normal',x);
    set(hpp(1),'Markersize',4,'marker','o','color','k')
    set(hpp(2),'Linestyle','-')
    annotation('textbox',[0.1582    0.4433    0.1066    0.0417],'string',['KS-p: ' num2str(KSp(1),3)],'Edgecolor','none');
    set(sp1,'Position',[0.1823    0.5838    0.2097    0.3412])
    %set(sp2,'Position',[0.4108    0.5838    0.2134    0.3412])
    set(sp3,'Position',[0.6502    0.5838    0.2097    0.3412])
    set(sp4,'Position',[0.1580    0.1387    0.2097    0.3412])
    if sum(x<=0) == 0
        sp5 = subplot(2,3,5);
        hpp = probplot('lognormal',x);
        set(hpp(1),'Markersize',4,'marker','o','color','k')
        set(hpp(2),'Linestyle','-')
        ylabel('')
        sp6 = subplot(2,3,6);
        hpp = probplot('weibull',x);
        set(hpp(1),'Markersize',4,'marker','o','color','k')
        set(hpp(2),'Linestyle','-')
        ylabel('')
        set(sp5,'Position',[0.4145    0.1400    0.2097    0.3412])
        set(sp6,'Position',[0.6693    0.1413    0.2097    0.3412])
        annotation('textbox',[0.4168    0.4407    0.0845    0.0417],'string',['KS-p: ' num2str(KSp(2),3)],'Edgecolor','none');
        annotation('textbox',[0.6651    0.4433    0.0953    0.0417],'string',['KS-p: ' num2str(KSp(3),3)],'Edgecolor','none');
    end
end