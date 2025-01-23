
%Example 3 of Cortez, Fluids 2021
%Channel with inflow and part of membrane permeable

%Developed by Ricardo Cortez, Brittany Leathers, and Michaela Kubacki
%July 2024

%Find best C1 for eps=Cds^(1/p) for different values of p
%And for different Flow Regimes

clear all 
% close all

%% Parameters to set 

%setting the viscosity
mu = 1; 
% Nvals=10*2.^(0:5);
Nvals=[10 20 40 80 160 320];
%number of points on boundary where velocity is set and force is computed 


%blob choice
blob=2;

%constant determining inflow velocity profile
a=4;

%channel
Lx = 5;
Ly = 1;
xmin = 0;
xmax = xmin + Lx;
ymin = 0;
ymax = ymin + Ly;
perm_min=5/3;  %start of permeable part
perm_max=10/3;  %end of permeable part


epsvals=0.00001*2.^(0:16);
epsvals=[epsvals (0.005:0.001:0.2)];
epsvals=sort(epsvals);
eumaxnorm=zeros(length(Nvals), length(epsvals));
eu2norm=zeros(length(Nvals), length(epsvals));

%%
for j=1:length(Nvals)
    N=Nvals(j)

%discretization of channel
ds = (ymax-ymin)/N;


for i=1:length(epsvals)
   

%%
%regularization parameter
ep=epsvals(i);

%discretization of top and bottom:
s = ds/2:ds:xmax-ds/2;
s = s';
%top wall coordinates (x,y) = (y1_top,y2_top)
y1_top=s;
y2_top = ymax*ones(size(y1_top));
%unit normals for top: 
normals_top=zeros(length(y1_top),2);
normals_top(:,2)=1;
%bottom wall coordinates (x,y) = (y1_bot,y2_bot)
y1_bot = s;
y2_bot = ymin*ones(size(y1_bot));
%unit normals for bottom
normals_bot=zeros(length(y1_bot),2); 
normals_bot(:,2)=-1;
%left-hand wall coordinates (x,y) = (y1_side,y2_side)
y2_side = (ymin+ds/2:ds:ymax-ds/2)';
y1_side = xmin*ones(size(y2_side));
%unit normals for side
normals_side=zeros(length(y1_side),2);
normals_side(:,1)=-1;

%coordiantes of boundary points
y1 = [y1_top;y1_bot; y1_side];
y2 = [y2_top;y2_bot; y2_side];
%For use at very end:
y2right = (ymin+ds/2:ds:ymax-ds/2)';
y1right = xmax*ones(size(y2_side));
y1f=[y1; y1right];
y2f=[y2; y2right];
%indices of permeable region
I=find(y1_top<2/3*Lx & y1_top>1/3*Lx);
I2=find(y1_top>=2/3*Lx | y1_top<=1/3*Lx);
%unit normals
normals=[normals_top; normals_bot; normals_side];




%velocity on boundary
%top
u1_top=zeros(size(y1_top));
u2_top=zeros(size(y1_top));
%bottom
u1_bot=zeros(size(y1_bot));
u2_bot=zeros(size(y1_bot));
%side: poiseulle flow
u1_side=a*(y2_side/Ly.*(1-y2_side/Ly));
u2_side=zeros(size(y1_side));

u1 = [u1_top; u1_bot; u1_side];
u2 = [u2_top; u2_bot; u2_side];


%computing the force 
f = RegStokeslets2D_velocitytoforce([y1,y2],[y1,y2],[u1,u2],ep,mu, blob);
f1 = f(:,1);
f2 = f(:,2);  

% Calculate velocities on right hand side too
Ufull=RegStokeslets2D_forcetovelocity([y1,y2],[f1,f2],[y1f,y2f],ep,mu, blob, ds);
ufull=Ufull(:,1); vfull=Ufull(:,2);




% Find solution at halfway mark
x2_half= (ymin+ds/2:ds:ymax-ds/2)';
x1_half= xmax/2*ones(size(x2_half));
uu = RegStokeslets2D_forcetovelocity([y1,y2],[f1,f2],...
    [x1_half,x2_half],ep,mu, blob);


uerror=abs(u1_side-uu(:,1));
eumaxnorm(j,i)=max(uerror);
eu2norm(j,i)=sqrt(ds*sum(uerror.^2));

end
end
%%
colorp=[0.4940, 0.1840, 0.5560];
colorlb=[0.3010, 0.6450, 0.9930];
colorg=[0.4660, 0.6740, 0.1880];
colordb=	[0, 0.4470, 0.7410];

forploty0=0.1;
forplotdx0=20;
dxforplot=epsvals;
yforplot=20*forploty0/forplotdx0^1*dxforplot.^1;
forploty02=20;
forplotdx02=10;
dxforplot2=epsvals;
yforplot2=forploty02/forplotdx02^2*dxforplot2.^(3/2);
figure;


   for j=1:length(Nvals)
   %Refinement study plots
   loglog(epsvals, abs(eumaxnorm(j,:)),'o-', 'LineWidth', 2.5, 'MarkerSize', ...
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
text(10^(-1),0.005,'$\epsilon$', 'Color', 'k', 'FontSize', 14, ...
    'Interpreter', 'latex');
hold off;
axis([10^(-6) 10^0 10^(-4) 10^(0) ]);
set(gca, 'FontSize', 17);
legend(sprintf('$N=%d$', Nvals(1)), sprintf('$N=%d$', Nvals(2)),...
    sprintf('$N=%d$', Nvals(3)), sprintf('$N=%d$', Nvals(4)),...
    sprintf('$N=%d$', Nvals(5)),sprintf('$N=%d$', Nvals(6)),...
    'Interpreter','latex', 'FontSize', 20);
xlabel('$\epsilon$', 'FontSize', 18,'Interpreter','latex')
ylabel('$||e||_{\infty}$', 'FontSize', 17,'Interpreter','latex')
title('Non-permeable channel: Error on line $x=2.5$',...
    'FontSize', 18,'Interpreter','latex')  
% yticks([10^(-5) 10^(-3) 10^(0)])
%%


colorp=[0.4940, 0.1840, 0.5560];
colorlb=[0.3010, 0.6450, 0.9930];
colorg=[0.4660, 0.6740, 0.1880];
colordb=	[0, 0.4470, 0.7410];

forploty0=0.1;
forplotdx0=20;
dxforplot=epsvals;
yforplot=20*forploty0/forplotdx0^1*dxforplot.^1;
forploty02=20;
forplotdx02=10;
dxforplot2=epsvals;
yforplot2=forploty02/forplotdx02^2*dxforplot2.^(3/2);
figure;
   for j=1:length(Nvals)
   %Refinement study plots
   
  
   loglog(epsvals, abs(eumaxnorm(j,:)),'o-', 'LineWidth', 2.5, 'MarkerSize', ...
       10)
hold on
% loglog(Nyvals, e2u, 's-', 'LineWidth', 2.5,'MarkerSize', 10, 'Color', colordb);
% loglog(Nyvals, e1u, 'd-', 'LineWidth', 2.5,'MarkerSize', 10,'Color', colorp);
% hold off
%legend('u', 'u away from boundary','u near boundary', 'FontSize', 12);
   end
   loglog(epsvals, yforplot, 'LineWidth', 2,'Color', 'k')
text(10^(-1),0.005,'$\epsilon$', 'Color', 'k', 'FontSize', 14, ...
    'Interpreter', 'latex');
hold off;
axis([10^(-6) 10^0 10^(-4) 10^(0) ]);
set(gca, 'FontSize', 17);
legend(sprintf('$N=%d$', Nvals(1)), sprintf('$N=%d$', Nvals(2)),...
    sprintf('$N=%d$', Nvals(3)), sprintf('$N=%d$', Nvals(4)),...
    sprintf('$N=%d$', Nvals(5)),sprintf('$N=%d$', Nvals(6)),...
    'Interpreter','latex', 'FontSize', 20);
xlabel('$\epsilon$', 'FontSize', 18,'Interpreter','latex')
ylabel('$||e||_{2}$', 'FontSize', 17,'Interpreter','latex')
title('Non-permeable channel: Error on line $x=2.5$',...
    'FontSize', 18,'Interpreter','latex')  



%%

besteps_max=zeros(length(Nvals),1);
besteps_2=zeros(length(Nvals),1);

for j=1:length(Nvals)

    [mine_max, kmax]=min(eumaxnorm(j,:));
    besteps_max(j)=epsvals(kmax);

    [mine_2, k2]=min(eu2norm(j,:));
    besteps_2(j)=epsvals(k2);

end


forploty0=0.1;
forplotdx0=1;
dxforplot=1./Nvals;
yforplot=forploty0/forplotdx0^1*dxforplot.^1;
figure;
loglog(1./Nvals, besteps_max,'o-', 'LineWidth', 2.5,...
    'MarkerSize', 10)
hold on;
loglog(1./Nvals, yforplot, 'LineWidth', 2,'Color', 'k')
text(10^(-2),0.002,'$\Delta s$', 'Color', 'k', 'FontSize', 14, ...
    'Interpreter', 'latex');
hold off;
axis([10^(-3) 0.2 10^(-5) 10^(0) ]);
set(gca, 'FontSize', 17);
xlabel('$\Delta s$', 'FontSize', 18,'Interpreter','latex')
ylabel('"Best" $\epsilon$ wrt max norm', 'FontSize', 17,'Interpreter','latex')
title('Non-permeable channel: ds vs epsilon',...
    'FontSize', 18,'Interpreter','latex')  



forploty0=0.1;
forplotdx0=1;
dxforplot=1./Nvals;
yforplot=forploty0/forplotdx0^1*dxforplot.^1;
figure;
loglog(1./Nvals, besteps_2,'o-', 'LineWidth', 2.5,...
    'MarkerSize', 10)
hold on;
loglog(1./Nvals, yforplot, 'LineWidth', 2,'Color', 'k')
text(10^(-2),0.002,'$\Delta s$', 'Color', 'k', 'FontSize', 14, ...
    'Interpreter', 'latex');
hold off;
axis([10^(-3) 0.2 10^(-5) 10^(0) ]);
set(gca, 'FontSize', 17);
xlabel('$\Delta s$', 'FontSize', 18,'Interpreter','latex')
ylabel('"Best" $\epsilon$ wrt 2 norm', 'FontSize', 17,'Interpreter','latex')
title('Non-permeable channel: ds vs epsilon',...
    'FontSize', 18,'Interpreter','latex')  


%% Best fit line
% Get coefficients of a line fit through the data.
coeffs_max = polyfit(1./Nvals, besteps_max, 1)
C0_max=coeffs_max(2);
C1_max=coeffs_max(1);
figure;
loglog(1./Nvals, besteps_max,'o-', 'LineWidth', 2.5,...
    'MarkerSize', 10)
hold on;
loglog(1./Nvals,C1_max./Nvals, 'LineWidth', 2,'Color', 'k')
hold off;
axis([10^(-3) 0.2 10^(-5) 10^(0) ]);
set(gca, 'FontSize', 17);
xlabel('$\Delta s$', 'FontSize', 18,'Interpreter','latex')
ylabel('"Best" $\epsilon$ wrt max norm', 'FontSize', 17,'Interpreter','latex')
title('Non-permeable channel: ds vs epsilon',...
    'FontSize', 18,'Interpreter','latex')  


% Get coefficients of a line fit through the data.
coeffs_2 = polyfit(1./Nvals, besteps_2, 1)
C0_2=coeffs_2(2);
C1_2=coeffs_2(1);

figure;
loglog(1./Nvals, besteps_2,'o-', 'LineWidth', 2.5,...
    'MarkerSize', 10)
hold on;
loglog(1./Nvals, C1_2./Nvals, 'LineWidth', 2,'Color', 'k')
hold off;
axis([10^(-3) 0.2 10^(-5) 10^(0) ]);
set(gca, 'FontSize', 17);
xlabel('$\Delta s$', 'FontSize', 18,'Interpreter','latex')
ylabel('"Best" $\epsilon$ wrt 2 norm', 'FontSize', 17,'Interpreter','latex')
title('Non-permeable channel: ds vs epsilon',...
    'FontSize', 18,'Interpreter','latex')  