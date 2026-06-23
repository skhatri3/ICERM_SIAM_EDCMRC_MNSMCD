clear all;

eps=[0.184 0.1 0.053 0.029 0.023 0.019 0.015 0.01 0.008 0.005];
dss=[1/10 1/20 1/40 1/80 1/100 1/120 1/160 1/240 1/320 1/500];

eps./(dss)
mean(eps(2:end)./(dss(2:end)))
eps./(dss.^(2/3))
mean(eps(4:end)./(dss(4:end).^(2/3)))
eps./(dss.^(1/2))
mean(eps(3:end)./(dss(3:end).^(1/2)))

%%

forploty0=0.3;
forplotdx0=5;
dxforplot=dss;
yforplot=20*forploty0/forplotdx0^1*dxforplot.^(2/3);

forploty0=0.3;
forplotdx0=5;
dxforplot=dss;
yforplot2=20*forploty0/forplotdx0^1*dxforplot.^(1);

forploty0=0.3;
forplotdx0=5;
dxforplot=dss;
yforplot3=20*forploty0/forplotdx0^1*dxforplot.^(1/2);

figure;
loglog(dss, eps,  'o-');
hold on;
loglog(dss, yforplot, 'LineWidth', 2,'Color', 'k')
loglog(dss, yforplot2, 'LineWidth', 2,'Color', 'k')
loglog(dss, yforplot3, 'LineWidth', 2,'Color', 'k')
text(0.01,0.01,'$\Delta s$', 'Color', 'k', 'FontSize', 14, ...
    'Interpreter', 'latex');
text(0.01,0.2,'$\Delta s^{2/3}$', 'Color', 'k', 'FontSize', 14, ...
    'Interpreter', 'latex');
text(0.006,0.06,'$\Delta s^{1/2}$', 'Color', 'k', 'FontSize', 14, ...
    'Interpreter', 'latex');
hold off;
% axis([10^(-6) 10^0 10^(-4) 10^(0) ]);
set(gca, 'FontSize', 17);
% legend(sprintf('$N=%d$', Nvals(1)), sprintf('$N=%d$', Nvals(2)),...
%     sprintf('$N=%d$', Nvals(3)), sprintf('$N=%d$', Nvals(4)),...
%     sprintf('$N=%d$', Nvals(5)),...
%     'Interpreter','latex', 'FontSize', 20);
xlabel('$\Delta s$', 'FontSize', 18,'Interpreter','latex')
ylabel('$\epsilon_{opt}$', 'FontSize', 17,'Interpreter','latex')
title('Non-permeable channel: On grid $[0.45, 0.55]\times [2,3]$',...
    'FontSize', 18,'Interpreter','latex')  

%%

p23=polyfit(dss.^(2/3), eps, 1)

p1=polyfit(dss, eps, 1)

p12=polyfit(dss.^(1/2), eps, 1)


%5.6 3/2
%%
p=polyfit(log(dss(7:10)), log(eps(7:10)), 1)

