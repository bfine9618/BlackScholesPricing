drift = log(close_price(2:end)) - log(close_price(1:end-1));
N = length(drift);
h = 1/365; %stock prices close one day apart

mu = sum(drift)/(N*h);
mu

var = sum((drift-mu*h).^2)/((N-1)*h);
var

%%
close all;
x = -.12:.01:.12;
n = histcounts(drift, x);

figure();
bar(x(1:end-1), n/N/.01);
xlim([-.12 .12]);
hold on;
norm = normpdf(x,mu*h,sqrt(var*h));
plot(x, norm, 'linewidth', 2, 'color','r')
title('Normalized Daily Drift vs. Expected Brownian');
saveas(gcf, 'brownianHist.png');

figure();
cdf = cumsum(n)/N;
bar(x(1:end-1), cdf);
hold on;
plot(x, normcdf(x, mu*h,sqrt(var*h)), 'linewidth', 2, 'color', 'r');
title('Normalized Daily Drift CDF vs. Expected Brownian CDF');
saveas(gcf, 'brownianCDF.png');

figure();
qqplot(drift);
title('QQ Plot to Test Normality');
saveas(gcf, 'qqplot.png');

%%

alph = .0375;

X0 = close_price(1);
Ex = X0*exp(mu+var/2);
K = [0.8 1 1.2]*Ex;
option_price = zeros(1, length(K));
risk_nuetral = alph - var/2;

for k=K
    a = (log(k/X0)-risk_nuetral)/sqrt(var);
    b = a - sqrt(var);
    QA = 1-normcdf(a, 0, 1);
    QB = 1 - normcdf(b, 0, 1);
    
    option_price(K==k) = X0*QB-exp(-alph)*k*QA;
end

disp(option_price);
