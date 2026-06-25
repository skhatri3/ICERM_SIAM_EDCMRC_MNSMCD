
%% Try one: using max norm, cutting off some points and looking at just 
% the part before hte minimum

clear all;
load C:\Users\britt\OneDrive\Documents\GitHub\ICERM_SIAM_EDCMRC_MNSMCD\Brittany\togeteps_opt.mat


eps=[0.184 0.1 0.053 0.029 0.023 0.015 0.01 0.008 0.005];



colorp=[0.4940, 0.1840, 0.5560];
colorlb=[0.3010, 0.6450, 0.9930];
colorg=[0.4660, 0.6740, 0.1880];
colordb=	[0, 0.4470, 0.7410];

forploty0=0.0001;
forplotdx0=1000;
dxforplot=epsvals;
yforplot=20*forploty0/forplotdx0^1.*dxforplot.^(2);
forploty02=20;
forplotdx02=10;
dxforplot2=epsvals;
yforplot2=forploty02/forplotdx02^2*dxforplot2.^(3/2);
figure;


for i=1:length(Nvals)

    slope = (log((eumaxnorm_away(i,2:end))./(eumaxnorm_away(i,1:end-1))))./...
        (log((epsvals(2:end))./(epsvals(1:end-1))));

   % opteps=eps(i);
   [~, Index]=find(slope<=-1 & slope>=-6 & epsvals(2:end)<=eps(i)  );
    Indices{i}=Index;
    slope(Index)
    %Refinement study plots
    loglog(epsvals(Index), eumaxnorm_away(i,Index),'o-', 'LineWidth', 2.5, 'MarkerSize', ...
        10)
    hold on;
    % loglog(Nyvals, e2u, 's-', 'LineWidth', 2.5,'MarkerSize', 10, 'Color', colordb);
    % loglog(Nyvals, e1u, 'd-', 'LineWidth', 2.5,'MarkerSize', 10,'Color', colorp);

    % loglog(epsvals, yforplot2, 'LineWidth', 2,'Color', 'b')
    % text(200,1.6/1000000,'$1/N^2$', 'Color', colorg, 'FontSize', 14, ...
    %     'Interpreter', 'latex');
    %legend('u', 'u away from boundary','u near boundary', 'FontSize', 12);
end
loglog(epsvals, yforplot, 'LineWidth', 2,'Color', 'k')
text(10^(-1),0.01,'$\epsilon$', 'Color', 'k', 'FontSize', 14, ...
    'Interpreter', 'latex');
hold off;
% axis([10^(-6) 10^0 10^(-4) 10^(0) ]);
set(gca, 'FontSize', 17);
legend(sprintf('$N=%d$', Nvals(1)), sprintf('$N=%d$', Nvals(2)),...
    sprintf('$N=%d$', Nvals(3)), sprintf('$N=%d$', Nvals(4)),...
    sprintf('$N=%d$', Nvals(5)),sprintf('$N=%d$', Nvals(6)),...
    sprintf('$N=%d$', Nvals(7)),sprintf('$N=%d$', Nvals(8)),...
    sprintf('$N=%d$', Nvals(9)),...
    'Interpreter','latex', 'FontSize', 20);
xlabel('$\epsilon$', 'FontSize', 18,'Interpreter','latex')
ylabel('$||e||_{\infty}$', 'FontSize', 17,'Interpreter','latex')
title('Non-permeable channel: Error on grid $[0.45, 0.55]\times [2,3]$',...
    'FontSize', 18,'Interpreter','latex')  
% yticks([10^(-5) 10^(-3) 10^(0)])
%$||e||_{\infty}$





NNs=9;
Neps=207;

% [Nvals, epsvals]=ndgrid(Nvals, epsvals);
% epsvals=reshape(epsvals, NNs*Neps, 1);
% Nvals=reshape(Nvals, NNs*Neps,1);
% 
% E=reshape(eumaxnorm_away, NNs*Neps,1);

epsvals_cut=[];
Nvals_cut=[];
E_cut=[];
for i=1:length(Nvals)
    epsvals_cut=[epsvals_cut; epsvals(Indices{i})'];
    Nvals_cut=[Nvals_cut; Nvals(i)*ones(length(Indices{i}), 1)];
    E_cut=[E_cut; eumaxnorm_away(i,Indices{i})'];
end


% Phi1 = [ (1./Nvals_cut(:)).^2./epsvals_cut(:).^2];
% Phi2 = [ (1./Nvals_cut(:)).^2./epsvals_cut(:)];
% Phi3 = [ (1./Nvals_cut(:)).^2];
% Phi4 = [ (1./Nvals_cut(:)).^2.*abs(log(epsvals_cut(:)))];

% Phi5= [  (1./Nvals_cut(:)).^6./epsvals_cut(:).^4];
% Phi6= [  (1./Nvals_cut(:)).^5./epsvals_cut(:).^3];
% Phi7= [  (1./Nvals_cut(:)).^4./epsvals_cut(:).^2];
% Phi8= [  (1./Nvals_cut(:)).^3./epsvals_cut(:).^1];

% coef1 = Phi2\E_cut(:)  
% res1 = norm(Phi1*coef1-E_cut(:))/norm(E_cut(:))
% coef2 = Phi2\E_cut(:)  
% res2 = norm(Phi2*coef2-E_cut(:))/norm(E_cut(:))
% coef3 = Phi3\E_cut(:)  
% res3 = norm(Phi3*coef3-E_cut(:))/norm(E_cut(:))
% coef4 = Phi4\E_cut(:) 
% res4 = norm(Phi4*coef4-E_cut(:))/norm(E_cut(:))

% coef5 = Phi5\E_cut(:) 
% res5 = norm(Phi5*coef5-E_cut(:))/norm(E_cut(:))
% coef6 = Phi6\E_cut(:) 
% res6 = norm(Phi6*coef6-E_cut(:))/norm(E_cut(:))
% coef7 = Phi7\E_cut(:) 
% res7 = norm(Phi7*coef7-E_cut(:))/norm(E_cut(:))
% coef8= Phi8\E_cut(:) 
% res8 = norm(Phi8*coef8-E_cut(:))/norm(E_cut(:))


%%  Try Two: Now use just one point so that the error is additive and only 
% cut off the beginning and try both errors together


clear all;
load C:\Users\britt\OneDrive\Documents\GitHub\ICERM_SIAM_EDCMRC_MNSMCD\Brittany\togeteps_opt.mat



eps_onepoint=[0.182 0.099 0.053 0.027 0.023 0.015 0.01 0.008 0.005];



colorp=[0.4940, 0.1840, 0.5560];
colorlb=[0.3010, 0.6450, 0.9930];
colorg=[0.4660, 0.6740, 0.1880];
colordb=	[0, 0.4470, 0.7410];

forploty0=0.0001;
forplotdx0=1000;
dxforplot=epsvals;
yforplot=20*forploty0/forplotdx0^1.*dxforplot.^(2);
forploty02=20;
forplotdx02=10;
dxforplot2=epsvals;
yforplot2=forploty02/forplotdx02^2*dxforplot2.^(3/2);
figure;


for i=1:length(Nvals)

    slope = (log((eu_pointaway(i,2:end))./(eu_pointaway(i,1:end-1))))./...
        (log((epsvals(2:end))./(epsvals(1:end-1))));

   % opteps=eps(i);
   [~, Index]=find(slope<=-0.025   );
    start{i}=min(Index);
   
    %Refinement study plots
    loglog(epsvals(start{i}:end), eumaxnorm_away(i,start{i}:end),'o-', 'LineWidth', 2.5, 'MarkerSize', ...
        10)
    hold on;
    % loglog(Nyvals, e2u, 's-', 'LineWidth', 2.5,'MarkerSize', 10, 'Color', colordb);
    % loglog(Nyvals, e1u, 'd-', 'LineWidth', 2.5,'MarkerSize', 10,'Color', colorp);

    % loglog(epsvals, yforplot2, 'LineWidth', 2,'Color', 'b')
    % text(200,1.6/1000000,'$1/N^2$', 'Color', colorg, 'FontSize', 14, ...
    %     'Interpreter', 'latex');
    %legend('u', 'u away from boundary','u near boundary', 'FontSize', 12);
end
loglog(epsvals, yforplot, 'LineWidth', 2,'Color', 'k')
text(10^(-1),0.01,'$\epsilon$', 'Color', 'k', 'FontSize', 14, ...
    'Interpreter', 'latex');
hold off;
% axis([10^(-6) 10^0 10^(-4) 10^(0) ]);
set(gca, 'FontSize', 17);
legend(sprintf('$N=%d$', Nvals(1)), sprintf('$N=%d$', Nvals(2)),...
    sprintf('$N=%d$', Nvals(3)), sprintf('$N=%d$', Nvals(4)),...
    sprintf('$N=%d$', Nvals(5)),sprintf('$N=%d$', Nvals(6)),...
    sprintf('$N=%d$', Nvals(7)),sprintf('$N=%d$', Nvals(8)),...
    sprintf('$N=%d$', Nvals(9)),...
    'Interpreter','latex', 'FontSize', 20);
xlabel('$\epsilon$', 'FontSize', 18,'Interpreter','latex')
ylabel('$||e||_{\infty}$', 'FontSize', 17,'Interpreter','latex')
title('Non-permeable channel: Error at one point away$',...
    'FontSize', 18,'Interpreter','latex')  
% yticks([10^(-5) 10^(-3) 10^(0)])
%$||e||_{\infty}$





NNs=9;
Neps=207;

% [Nvals, epsvals]=ndgrid(Nvals, epsvals);
% epsvals=reshape(epsvals, NNs*Neps, 1);
% Nvals=reshape(Nvals, NNs*Neps,1);
% 
% E=reshape(eumaxnorm_away, NNs*Neps,1);

epsvals_cut= [];
Nvals_cut=[];
E_cut=[];
for i=1:length(Nvals)
    epsvals_cut=[epsvals_cut; epsvals(start{i}:end)'];
    Nvals_cut=[Nvals_cut; Nvals(i)*ones(length(epsvals(start{i}:end)), 1)];
    E_cut=[E_cut; eu_pointaway(i,start{i}:end)'];
end

% 
% Phi1 = [epsvals_cut(:).^2, (1./Nvals_cut(:)).^2./epsvals_cut(:).^2];
% Phi2 = [epsvals_cut(:).^2, (1./Nvals_cut(:)).^2./epsvals_cut(:)];
% Phi3 = [epsvals_cut(:).^2, (1./Nvals_cut(:)).^2];
% Phi4 = [epsvals_cut(:).^2, (1./Nvals_cut(:)).^2.*abs(log(epsvals_cut(:)))];
% 
% coef1 = Phi1\E_cut(:) 
% res1 = norm(Phi1*coef1-E_cut(:))/norm(E_cut(:))
% coef2 = Phi2\E_cut(:)  
% res2 = norm(Phi2*coef2-E_cut(:))/norm(E_cut(:))
% coef3 = Phi3\E_cut(:)  
% res3 = norm(Phi3*coef3-E_cut(:))/norm(E_cut(:))
% coef4 = Phi4\E_cut(:) 
% res4 = norm(Phi4*coef4-E_cut(:))/norm(E_cut(:))

Phi5= [  (1./Nvals_cut(:)).^6./epsvals_cut(:).^4];
Phi6= [  (1./Nvals_cut(:)).^5./epsvals_cut(:).^3];
Phi7= [  (1./Nvals_cut(:)).^4./epsvals_cut(:).^2];
Phi8= [  (1./Nvals_cut(:)).^3./epsvals_cut(:).^1];

coef5 = Phi5\E_cut(:) 
res5 = norm(Phi5*coef5-E_cut(:))/norm(E_cut(:))
coef6 = Phi6\E_cut(:) 
res6 = norm(Phi6*coef6-E_cut(:))/norm(E_cut(:))
coef7 = Phi7\E_cut(:) 
res7 = norm(Phi7*coef7-E_cut(:))/norm(E_cut(:))
coef8= Phi8\E_cut(:) 
res8 = norm(Phi8*coef8-E_cut(:))/norm(E_cut(:))


% Phi=[epsvals_cut(:).^2, (1./Nvals_cut(:)).^2./epsvals_cut(:).^2,...
%     (1./Nvals_cut(:)).^2./epsvals_cut(:), (1./Nvals_cut(:)).^2, ...
%     (1./Nvals_cut(:)).^2.*abs(log(epsvals_cut(:)))];
% 
% coef = Phi\E_cut(:) 
% res = norm(Phi*coef-E_cut(:))/norm(E_cut(:))

%% Try three: Just look at a plot of the slopes before the minimum and see 
%we can nail down k_3


clear all;
load C:\Users\britt\OneDrive\Documents\GitHub\ICERM_SIAM_EDCMRC_MNSMCD\Brittany\togeteps_opt.mat

eps=[0.184 0.1 0.053 0.029 0.023 0.015 0.01 0.008 0.005];
eps_onepoint=[0.182 0.099 0.053 0.027 0.023 0.015 0.01 0.008 0.005];

epsvals=unique(epsvals);

for i=1:length(Nvals)

    opteps=eps(i);
    Ind=find(abs(epsvals-opteps)<10^(-6));

    % epsvals(Ind)

    % figure;
    % loglog(epsvals(1:Ind), eumaxnorm_away(i,1:Ind),'o-', 'LineWidth', 2.5, 'MarkerSize', ...
    %     10)

    slope = (log((eumaxnorm_away(i,2:Ind))./(eumaxnorm_away(i,1:Ind-1))))./...
        (log((epsvals(2:Ind))./(epsvals(1:Ind-1))));

    figure; plot(epsvals(2:Ind), slope, 'o');
    title('Fixed N, slope of loglog plot of epsilon vs. Max Error')
    xlabel('\epsilon', 'Interpreter','latex')
    ylabel('$||e||_{\infty}$', 'Interpreter','latex')



    % figure;
    % edges = -5:0.5:5;
    % histogram(slope,edges)
end