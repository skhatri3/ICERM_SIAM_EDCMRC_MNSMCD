%Example 3 of Cortez, Fluids 2021
%Channel with inflow and part of membrane permeable

%Developed by Ricardo Cortez, Brittany Leathers, and Michaela Kubacki
%July 2024


clear all 
% close all

%% Parameters to set 

%setting the viscosity
mu = 1; 

%number of points on boundary where velocity is set and force is computed 
N = 160;    

%blob choice
blob=2;

%permeability coefficient (assuming constant for the permeable region)
b=-0.00011;

%constant determining inflow velocity profile
a=4;

%% Setting forces and computing velocity 
%channel
Lx = 5;
Ly = 1;
xmin = 0;
xmax = xmin + Lx;
ymin = 0;
ymax = ymin + Ly;
perm_min=5/3;  %start of permeable part
perm_max=10/3;  %end of permeable part

%discretization of channel
ds = (ymax-ymin)/N;
%%
%regularization parameter
ep =0.0224;%1*ds;

%discretization of top:
s1=(xmin:ds:perm_min-ds)';
s_perm=(perm_min:ds:perm_max)';
s2=(perm_max+ds:ds:xmax)';
%top wall coordinates (x,y) = (y1_top,y2_top)
y1_top1=s1;
y2_top1 = ymax*ones(size(y1_top1));
y1_top_perm=s_perm;
y2_top_perm = ymax*ones(size(y1_top_perm));
y1_top2=s2;
y2_top2 = ymax*ones(size(y1_top2));

%unit normals for top: 
normals_top=zeros(length([y1_top1; y1_top_perm; y1_top2]),2);
normals_top(:,2)=1;

%discretization of bottom and left:
s = 0:ds:xmax;
s = s';
%bottom wall coordinates (x,y) = (y1_bot,y2_bot)
y1_bot = s;
y2_bot = ymin*ones(size(y1_bot));
%unit normals for bottom
normals_bot=zeros(length(y1_bot),2); 
normals_bot(:,2)=-1;
%left-hand wall coordinates (x,y) = (y1_side,y2_side)
y2_side = (ymin+ds:ds:ymax-ds)';
y1_side = xmin*ones(size(y2_side));
%unit normals for side
normals_side=zeros(length(y1_side),2);
normals_side(:,1)=-1;

%coordiantes of boundary points
y1 = [y1_top1; y1_top_perm;y1_top2;y1_bot; y1_side];
y2 = [y2_top1; y2_top_perm;y2_top2;y2_bot; y2_side];
%indices of permeable region
k1=length(y1_top1);
k2=length(y1_top_perm);
I=((k1+1):(k1+k2))';
%unit normals
normals=[normals_top; normals_bot; normals_side];
%beta vector for function;
beta= zeros(length(y1),1); 
beta(I)=b;



%velocity on boundary
%top
u1_top1=zeros(size(y1_top1));
u2_top1=zeros(size(y1_top1));
u1_top2=zeros(size(y1_top2));
u2_top2=zeros(size(y1_top2));
%permeable part : temporary velocity for finding g
u1_top_perm=zeros(size(y1_top_perm));
u2_top_perm=zeros(size(y1_top_perm));
%bottom
u1_bot=zeros(size(y1_bot));
u2_bot=zeros(size(y1_bot));
%side: poiseulle flow
u1_side=a*(y2_side/Ly.*(1-y2_side/Ly));
u2_side=zeros(size(y1_side));

u1 = [u1_top1;u1_top_perm;u1_top2; u1_bot; u1_side];
u2 = [u2_top1;u2_top_perm;u2_top2; u2_bot; u2_side];



%compute g
g=RegStokeslets2D_velocityto_gforce_permeable([y1,y2],[y1,y2],...
    [u1,u2],ep,mu, blob, I, beta, normals);

%Find velocity in permeable region:
y1b=y1(I);
y2b=y2(I);
[u_beta]=RegStokeslets2D_permeable_gtovelocity([y1,y2],g, [y1b,y2b],...
    ep,mu, blob, beta, normals);

%Put in the new velocities for the permeable part
%permeable part : temporary velocity for finding g
u1(I)=u_beta(:,1);
u2(I)=u_beta(:,2);

%computing the force 
f = RegStokeslets2D_velocitytoforce([y1,y2],[y1,y2],[u1,u2],ep,mu, blob);
f1 = f(:,1);
f2 = f(:,2);  




%%

%% Compute velocity everywhere on a grid to view
%domain on which velocity is computed
x1min = 0;
x1max = xmax;
x2min = 0;
x2max = ymax;
%resolution for velocity grid
Nx1 = 100;
Nx2 = 100;
%points on which velocity will be computed
xx1 = linspace(x1min,x1max,Nx1);
xx2 = linspace(x2min,x2max,Nx2);
[x1m,x2m] = meshgrid(xx1,xx2);
x1 = x1m(:);
x2 = x2m(:);
%computing velocity
uu = RegStokeslets2D_forcetovelocity([y1,y2],[f1,f2],[x1,x2],ep,mu, blob);
uu1 = uu(:,1);
uu2 = uu(:,2);
u1m = reshape(uu1,size(xx1,2),size(xx2,2));
u2m = reshape(uu2,size(xx1,2),size(xx2,2));
% Plot channel with velocity everywhere in channel
figure;
plot(y1,y2,'k.')
pcolor(x1m, x2m, u1m)
% box on
hold on
    h = streamslice(x1m,x2m,u1m,u2m);
    set( h, 'Color','k' )
    set( h, 'LineWidth', 1)
colorbar
shading flat
axis equal
xlim([xmin,xmax])
ylim([ymin,ymax])
title('Velocity Magnitude')
hold off



% %% Plotting figures 
% set(0,'defaultaxesfontsize',20,'defaultaxeslinewidth',2.0,...
%       'defaultlinelinewidth',2.0,'defaultlinemarkersize',10.0)
% 
% figure;
% plot(y1,y2, 'o');
% hold on;
% quiver(y1,y2, f1, f2)
% axis([xmin-0.5 xmax+0.5 ymin-0.5 ymax+0.5])


%% Check flow rates 

% inlet
Rin=ds*sum(dot(normals_side, [u1_side u2_side]))

%top
normals_top_perm=zeros(length( y1_top_perm),2);
normals_top_perm(:,2)=1;
Rtop=ds*sum(dot(normals_top_perm, u_beta))

%outlet
x2_out= (ymin+ds:ds:ymax-ds)';
x1_out= xmax*ones(size(x2_out));
uu = RegStokeslets2D_forcetovelocity([y1,y2],[f1,f2],...
    [x1_out,x2_out],ep,mu, blob);
Rout=-ds*sum(dot(normals_side, uu))

Rin+Rtop+Rout