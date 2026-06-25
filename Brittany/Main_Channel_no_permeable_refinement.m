
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
Nvals=[10 20 40 80 160 320 640 ];
%number of points on boundary where velocity is set and force is computed 


%blob choice
blob=2;

%constant determining inflow velocity profile
a=4;
%% 

%channel
Lx = 5;
Ly = 1;
xmin = 0;
xmax = xmin + Lx;
ymin = 0;
ymax = ymin + Ly;

eumaxnorm_away=zeros(length(Nvals),1);
eu2norm_away=zeros(length(Nvals), 1);
eu_pointaway=zeros(length(Nvals), 1);
eumaxnorm_near=zeros(length(Nvals),1);
eu2norm_near=zeros(length(Nvals), 1);

%% Refinement Study

for i=1:length(Nvals)
    N=Nvals(i)

%discretization of channel
ds = (ymax-ymin)/N;
ds_x=ds; ds_y=ds;

%Regularization parameter
% ep=1.6645*ds;
% ep=1.8*ds;
% 
% % ep=2*sqrt(ds)
% ep=0.9*ds^(2/3);
% ep=0.35*ds^(2/3);
ep=0.15*ds^(1/2);

% ep=2.5*ds;

ep=2.7*ds;
% ep=0.43*ds^(2/3);
ep=0.2*ds^(1/2);
ep=0.1*ds^(1/2);


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

% Calculate velocities on right hand side too
Ufull=RegStokeslets2D_forcetovelocity([y1,y2],[f1,f2],[y1f,y2f],ep,mu,...
    blob, wt);
ufull=Ufull(:,1); vfull=Ufull(:,2);


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
eumaxnorm_away(i)=max(uerror_away(:));
eu2norm_away(i)=sqrt(dx_away^2*sum(sum(uerror_away.^2)));

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
eumaxnorm_near(i)=max(uerror_near(:));
eu2norm_near(i)=sqrt(dx_near^2*sum(sum(uerror_near.^2)));



%Find solution at single point away from boundary and single point near
%boundary
%Away: (2,1/2)
uaway=RegStokeslets2D_forcetovelocity([y1,y2],[f1,f2],[2 0.5],ep,mu,...
    blob, wt);
eu_pointaway(i)=abs(uaway(:,1)-a*(0.5/Ly.*(1-0.5/Ly)));

% % Calculate velocities on right hand side too
% Ufull=RegStokeslets2D_forcetovelocity([y1,y2],[f1,f2],[y1f,y2f],ep,mu,...
%     blob, wt);
% ufull=Ufull(:,1); vfull=Ufull(:,2);
% 
% % Check flow rates 
% % inlet
% Rin=ds*sum(dot(normals_side, [u1_side u2_side]));
% 
% %top
% normals_top=zeros(length( y1_top),2);
% normals_top(:,2)=1;
% Rtop=ds*sum(dot(normals_top, [ufull(1:length(y1_top)), ...
%     vfull(1:length(y1_top))]));
% 
% %outlet
% x2_out= (ymin+ds/2:ds:ymax-ds/2)';
% x1_out= xmax*ones(size(x2_out));
% uu = RegStokeslets2D_forcetovelocity([y1,y2],[f1,f2],...
%     [x1_out,x2_out],ep,mu, blob, wt);
% Rout=-ds*sum(dot(normals_side, uu));
% 
% error(i)=Rin+Rtop+Rout;


end
%%

colorp=[0.4940, 0.1840, 0.5560];
colorlb=[0.3010, 0.6450, 0.9930];
colorg=[0.4660, 0.6740, 0.1880];
colordb=	[0, 0.4470, 0.7410];

forploty0=0.1;
forplotdx0=0.02;
dxforplot=1./Nvals;
yforplot=forploty0/forplotdx0^1*dxforplot.^(2);
forploty02=25;
forplotdx02=10;
dxforplot2=1./Nvals;
yforplot2=forploty02/forplotdx02^2*dxforplot2.^1;
   
%Refinement study plots
figure;
loglog(Nvals, eumaxnorm_near,'o-', 'LineWidth', 2.5, 'MarkerSize', ...
       10, 'Color', colorlb)
hold on
% loglog(Nyvals, e2u, 's-', 'LineWidth', 2.5,'MarkerSize', 10, 'Color', colordb);
% loglog(Nyvals, e1u, 'd-', 'LineWidth', 2.5,'MarkerSize', 10,'Color', colorp);
   loglog(Nvals, yforplot, 'LineWidth', 2,'Color', colorg)
text(100,1/7000,'$\Delta s^{2}$', 'Color', colorg, 'FontSize', 14, ...
    'Interpreter', 'latex');
% loglog(Nvals, yforplot2, 'LineWidth', 2,'Color', colorg)
% text(100,1/300,'$\Delta s$', 'Color', colorg, 'FontSize', 14, ...
%     'Interpreter', 'latex');
hold off
%legend('u', 'u away from boundary','u near boundary', 'FontSize', 12);
axis([10 1000 10^(-7) 10^(0) ]);
set(gca, 'FontSize', 17);
% legend('$L^{\infty}$', '$L^2$', '$L^1$','Interpreter','latex', 'FontSize', 20);
xlabel('$N$', 'FontSize', 18,'Interpreter','latex')
ylabel('$||e||_{\infty}$', 'FontSize', 17, 'Interpreter','latex')
title('Refinement Study, $\epsilon=0.1ds^{1/2}$, "away"',...
    'FontSize', 18,'Interpreter','latex')  
% yticks([10^(-5) 10^(-3) 10^(0)])


% %%
% d1u=zeros(1, length(Nvals));
% d2u=zeros(1, length(Nvals));
% dmaxu=zeros(1, length(Nvals));
% 
% d1v=zeros(1, length(Nvals));
% d2v=zeros(1, length(Nvals));
% dmaxv=zeros(1, length(Nvals));
% 
% for i=1:length(Nvals)-1
%     %from fine to coarse 
%     ufine=usolutions{i+1}; 
%     ucoarse=usolutions{i};
%     diffu=(ucoarse-ufine);
%     % ufinerest=ufine(1:2:end, 1:2:end);
%     % diffu=(ucoarse-ufinerest);
% 
% 
%     vfine=vsolutions{i+1}; 
%     vcoarse=vsolutions{i};
%     diffv=(vcoarse-vfine);
%     % vfinerest=vfine(1:2:end, 1:2:end);
%     % diffv=(vcoarse-vfinerest);
%     % dx=dxs(i);
% 
%     d1u(1, i)= dx^2*sum(sum(abs(diffu)));
%     d2u(1, i)= sqrt(dx^2*sum(sum(diffu(:,:).^2)));
%     dmaxu(1, i)= max(max(abs(diffu)));
% 
%     d1v(1, i)= dx^2*sum(sum(abs(diffv)));
%     d2v(1, i)= sqrt(dx^2*sum(sum(diffv(:,:).^2)));
%     dmaxv(1, i)= max(max(abs(diffv)));
% 
% end

%%

% %% Allnorms refinement study
% 
% colorp=[0.4940, 0.1840, 0.5560];
% colorlb=[0.3010, 0.6450, 0.9930];
% colorg=[0.4660, 0.6740, 0.1880];
% colordb=	[0, 0.4470, 0.7410];
% 
% forploty0=1.5;
% forplotdx0=10;
% dxforplot=1./Nvals;
% yforplot=forploty0/forplotdx0^1*dxforplot.^1;
% forploty02=50;
% forplotdx02=10;
% dxforplot2=1./Nvals;
% yforplot2=forploty02/forplotdx02^2*dxforplot2.^2;
% 
% %Refinement study plots
% figure;
% loglog(Nvals, dmaxu,'o-', 'LineWidth', 3, 'MarkerSize', ...
%        10, 'Color', colorlb)
% hold on
% loglog(Nvals, d2u, 's-', 'LineWidth', 3,'MarkerSize', 10, 'Color', colordb);
% loglog(Nvals, d1u, 'd-', 'LineWidth', 3,'MarkerSize', 10,'Color', colorp);
%    loglog(Nvals(3:5), yforplot(3:5), 'LineWidth', 2,'Color', colorg)
% text(360,1/7000,'1/N', 'Color', colorg, 'FontSize', 16);
% % loglog(Nvals(3:5), yforplot2(3:5), 'LineWidth', 1.5,'Color', colorg)
% % text(100,1.4/100000,'1/N^2', 'Color', colorg, 'FontSize', 12);
% hold off
% %legend('u', 'u away from boundary','u near boundary', 'FontSize', 12);
% legend('L^{\infty}', 'L^2', 'L^1');
% axis([10 360 10^(-4) 10^(0)]);
% set(gca, 'FontSize', 12);
% xlabel('N_x', 'FontSize', 18)
% ylabel('Difference Norms', 'FontSize', 18)
% title('Horizontal Error Refinement Study',...
%     'FontSize', 18,'Interpreter','latex')
   