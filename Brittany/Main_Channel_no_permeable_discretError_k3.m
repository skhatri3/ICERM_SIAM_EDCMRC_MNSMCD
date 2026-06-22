
%Example 3 of Cortez, Fluids 2021
%Channel with inflow and part of membrane permeable

%Developed by Ricardo Cortez, Brittany Leathers, and Michaela Kubacki
%July 2024


clear all 
% close all

%% Parameters to set 

%setting the viscosity
mu = 1; 
% Nvals=10*2.^(0:5);
Nvals=[10 20 40 80 100 160 240 320 500];
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

%Regularization parameter
% epsvals=0.00001*2.^(0:16);
% epsvals=[epsvals (0.05:0.03:0.2) 0.25 0.4 0.5 0.6];
% epsvals=sort(epsvals);

epsvals=[0.003:0.001:0.015, 0.015:0.002:0.1 0.1:0.002:0.4];

% nearpointoptions=[0.001 0.005 0.01 0.05 0.1];

eumaxnorm_away=zeros(length(Nvals), length(epsvals));
eu2norm_away=zeros(length(Nvals), length(epsvals));
eumaxnorm_near=zeros(length(Nvals), length(epsvals));
eu2norm_near=zeros(length(Nvals), length(epsvals));
eu_pointaway=zeros(length(Nvals), length(epsvals));
% eu_pointnear=zeros(length(Nvals), length(epsvals), length(nearpointoptions));

%% Loop through Nvals

for i=1:length(Nvals)
    N=Nvals(i)



%discretization of channel
ds = (ymax-ymin)/N;
ds_x=ds; ds_y=ds;

%% Loop through epsilon values

for j=1:length(epsvals)

    %Regularization parameter
    ep=epsvals(j);
   

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

% % %unit normals
% % normals=[normals_top; normals_bot; normals_side];

%weights
wt = [ds_x*ones(size(y1_top)); ds_x*ones(size(y1_bot));...
    ds_y*ones(size(y1_side))];



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
f = RegStokeslets2D_velocitytoforce([y1,y2],[y1,y2],[u1,u2],ep,mu, blob,...
    wt);
f1 = f(:,1);
f2 = f(:,2);  


%calculate on grid away from the boundary
xmin_away=2;
ymin_away=0.45;
Lx_away=1; xmax_away=xmin_away+Lx_away;
Ly_away=0.1; ymax_away=ymin_away+Ly_away;
dx_away=(ymax_away-ymin_away)/(2^6);
% dx=ds;
% dxs(i)=dx;
Nx_away=round(Lx_away/dx_away);
Ny_away=round(Ly_away/dx_away);
xg_away=dx_away*(0:Nx_away-1)+xmin_away;
yg_away=dx_away*(0:Ny_away-1)+ymin_away;
[xg_away,yg_away]=ndgrid(xg_away,yg_away);
xgv_away=reshape(xg_away, Nx_away*Ny_away,1);
ygv_away=reshape(yg_away, Nx_away*Ny_away,1);
Ugrid=RegStokeslets2D_forcetovelocity([y1,y2],[f1,f2],[xgv_away,ygv_away],ep,mu,...
    blob, wt);
ug_away=reshape(Ugrid(:,1), Nx_away, Ny_away); 
vg_away=reshape(Ugrid(:,2), Nx_away, Ny_away); 

u_soln_away=a*(yg_away/Ly.*(1-yg_away/Ly));
uerror_away=abs(ug_away-u_soln_away);
eumaxnorm_away(i,j)=max(uerror_away(:));
eu2norm_away(i,j)=sqrt(dx_away^2*sum(sum(uerror_away.^2)));
% 
%calculate on grid including near the boundary
xmin_near=2;
ymin_near=0.01;
Lx_near=1; xmax_near=xmin_near+Lx_near;
Ly_near=0.4; ymax_near=ymin_near+Ly_near;
dx_near=(ymax_near-ymin_near)/(2^6);
% dx=ds;
% dxs(i)=dx;
Nx_near=round(Lx_near/dx_near);
Ny_near=round(Ly_near/dx_near);
xg_near=dx_near*(0:Nx_near-1)+xmin_near;
yg_near=dx_near*(0:Ny_near-1)+ymin_near;
[xg_near,yg_near]=ndgrid(xg_near,yg_near);
xgv_near=reshape(xg_near, Nx_near*Ny_near,1);
ygv_near=reshape(yg_near, Nx_near*Ny_near,1);
Ugrid=RegStokeslets2D_forcetovelocity([y1,y2],[f1,f2],[xgv_near,ygv_near],ep,mu,...
    blob, wt);
ug_near=reshape(Ugrid(:,1), Nx_near, Ny_near); 
vg_near=reshape(Ugrid(:,2), Nx_near, Ny_near); 

u_soln_near=a*(yg_near/Ly.*(1-yg_near/Ly));
uerror_near=abs(ug_near-u_soln_near);
eumaxnorm_near(i,j)=max(uerror_near(:));
eu2norm_near(i,j)=sqrt(dx_near^2*sum(sum(uerror_near.^2)));


%Find solution at single point away from boundary and single point near
%boundary
%Away: (2,1/2)
uaway=RegStokeslets2D_forcetovelocity([y1,y2],[f1,f2],[2 0.5],ep,mu,...
    blob, wt);
eu_pointaway(i,j)=abs(uaway(:,1)-a*(0.5/Ly.*(1-0.5/Ly)));

% %Near: 
% for k=1:length(nearpointoptions)
%     if k==1
%         unear=RegStokeslets2D_forcetovelocity([y1,y2],[f1,f2],[2 sqrt(ep)/2],...
%             ep,mu,blob, wt);
%         eu_pointnear(i,j, k)=abs(unear(:,1)-a*(sqrt(ep)/2/Ly.*(1-sqrt(ep)/2/Ly)));
%     else
%         unear=RegStokeslets2D_forcetovelocity([y1,y2],[f1,f2],...
%             [2 nearpointoptions(k)],...
%             ep,mu,blob, wt);
%         eu_pointnear(i,j, k)=abs(unear(:,1)-a*(nearpointoptions(k)/Ly.*...
%             (1-nearpointoptions(k)/Ly)));
%     end
% end

end

end

%%

% clear all
% 
% Data=open('Error_eps_no_perm.mat');
% Nvals=Data.Nvals;
% epsvals=Data.epsvals;
% eu2norm_away=Data.eu2norm_away;
% eu2norm_near=Data.eu2norm_near;
% eu_pointaway=Data.eu_pointaway;
% eumaxnorm_away=Data.eumaxnorm_away;
% eumaxnorm_near=Data.eumaxnorm_near;
% 
% clear('Data');

%% Plot error

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

   %Refinement study plots
   loglog(epsvals, eu_pointaway(i,:),'o-', 'LineWidth', 2.5, 'MarkerSize', ...
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



