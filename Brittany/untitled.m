%% Plot error

colorp=[0.4940, 0.1840, 0.5560];
colorlb=[0.3010, 0.6450, 0.9930];
colorg=[0.4660, 0.6740, 0.1880];
colordb=	[0, 0.4470, 0.7410];

forploty0=0.1;
forplotdx0=15;
dxforplot=epsvals;
yforplot=20*forploty0/forplotdx0^1*dxforplot.^(2);
forploty02=20;
forplotdx02=10;
dxforplot2=epsvals;
yforplot2=forploty02/forplotdx02^2*dxforplot2.^(3/2);
figure;


   for i=1:length(Nvals)

   %Refinement study plots
   loglog(epsvals, eu_pointnear(i,:, 3),'o-', 'LineWidth', 2.5, 'MarkerSize', ...
       10)
    hold on
% loglog(Nyvals, e2u, 's-', 'LineWidth', 2.5,'MarkerSize', 10, 'Color', colordb);
% loglog(Nyvals, e1u, 'd-', 'LineWidth', 2.5,'MarkerSize', 10,'Color', colorp);

% loglog(epsvals, yforplot2, 'LineWidth', 2,'Color', 'b')
% text(200,1.6/1000000,'$1/N^2$', 'Color', colorg, 'FontSize', 14, ...
%     'Interpreter', 'latex');
%legend('u', 'u away from boundary','u near boundary', 'FontSize', 12);
   end
loglog(epsvals, yforplot, 'LineWidth', 2,'Color', 'k')
text(10^(-2),0.00001,'$\epsilon^2$', 'Color', 'k', 'FontSize', 14, ...
    'Interpreter', 'latex');
hold off;
axis([10^(-6) 10^0 10^(-4) 10^(0) ]);
set(gca, 'FontSize', 17);
legend(sprintf('$N=%d$', Nvals(1)), sprintf('$N=%d$', Nvals(2)),...
    'Interpreter','latex', 'FontSize', 20);
xlabel('$\epsilon$', 'FontSize', 18,'Interpreter','latex')
ylabel('Error', 'FontSize', 17,'Interpreter','latex')
title('Non-permeable channel: Error on point $(2,0.01)$',...
    'FontSize', 18,'Interpreter','latex')  
% yticks([10^(-5) 10^(-3) 10^(0)])
%$||e||_{\infty}$
%sprintf('$N=%d$', Nvals(1)), sprintf('$N=%d$', Nvals(2)),...
%sprintf('$N=%d$', Nvals(3)), sprintf('$N=%d$', Nvals(4)),...
 %   sprintf('$N=%d$', Nvals(5)),...


