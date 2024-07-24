

%Poiseulle flow 

%Developed by Shilpa Khatri and Ricardo Cortez  
%July 2024 

clear all 
% close all

%% Parameters to set 

%setting the viscosity
mu = 1; 

%number of points on boundary where velocity is set and force is computed 
N = 600;    

%blob choice
blob=3;

%% Setting forces and computing velocity 



%channel
Lx = 8;
Ly = 1;
xmin = 0;
xmax = xmin + Lx;
ymin = 0;
ymax = ymin + Ly;

%discretization of channel
ds = xmax/N;
s = 0:ds:xmax;
s = s';
%regularization parameter
ep = 1.75*ds;
%discretization of top:
s1=(0:ds:1-ds)';
s2=(1:ds:2)';
s3=(2+ds:ds:xmax)';

%top wall coordinates (x,y) = (y1_top,y2_top)
y1_top1=s1;
y2_top1 = ymax*ones(size(y1_top1));
y1_top2=s2;
y2_top2 = ymax*ones(size(y1_top2));
y1_top3=s3;
y2_top3 = ymax*ones(size(y1_top3));


%bottom wall coordinates (x,y) = (y1_bot,y2_bot)
y1_bot = s;
y2_bot = ymin*ones(size(y1_bot));
%left-hand wall coordinates (x,y) = (y1_side,y2_side)
y2_side = (ymin+ds:ds:ymax-ds)';
y1_side = zeros(size(y2_side));
%coordiantes of boundary points
y1 = [y1_top1; y1_top2;y1_top3;y1_bot; y1_side];
y2 = [y2_top1; y2_top2;y2_top3;y2_bot; y2_side];



%velocity on boundary
%top
u1_top1=zeros(size(y1_top1));
u2_top1=zeros(size(y1_top1));
u1_top3=zeros(size(y1_top3));
u2_top3=zeros(size(y1_top3));
%permeable part 
%%%% THIS PART WON"T BE TRUE ONCE YOU PUT IN PERMEABLILITY
u1_top2=zeros(size(y1_top2));
u2_top2=zeros(size(y1_top2));
%bottom
u1_bot=zeros(size(y1_bot));
u2_bot=zeros(size(y1_bot));
%side
a=10;
u1_side=-a*(y2_side-ymin).*(y2_side-ymax);
u2_side=zeros(size(y1_side));

u1 = [u1_top1;u1_top2;u1_top3; u1_bot; u1_side];
u2 = [u2_top1;u2_top2;u2_top3; u2_bot; u2_side];


%computing the force 
f = RegStokeslets2D_velocitytoforce([y1,y2],[y1,y2],[u1,u2],ep,mu, blob);

%f is a force and to compare to exact solution need force density - divide
%by radius*dt 
f1 = f(:,1);
f2 = f(:,2);  

%find velocities to the right in the channel
x2=y2_side;
x1=5*ones(size(x2));
u = RegStokeslets2D_forcetovelocity([y1,y2],f,[x1,x2],ep,mu,blob);

error1=abs(u(:,1)-u1_side);
error2=abs(u(:,2)-u2_side);

max(max(error1))/max(u1_side)
max(max(error2))

figure; plot(x1, x2); hold on; quiver(u(:,1), u(:,2))



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
[x1m,x2m] = ndgrid(xx1,xx2);
x1 = x1m(:);
x2 = x2m(:);
%computing velocity
uu = RegStokeslets2D_forcetovelocity([y1,y2],[f1,f2],[x1,x2],ep,mu, blob);
uu1 = uu(:,1);
uu2 = uu(:,2);
u1m = reshape(uu1,size(xx1,2),size(xx2,2));
u2m = reshape(uu2,size(xx1,2),size(xx2,2));
% Plot channel with velocity everywhere in channel
figure(2)
plot(y1,y2,'k.')
surf(x1m, x2m, u1m, 'EdgeColor', 'none')
box on
colorbar
view(2)
axis equal
xlim([xmin,xmax])
ylim([ymin,ymax])
title('Velocity Magnitude')




%% Plotting figures 
set(0,'defaultaxesfontsize',20,'defaultaxeslinewidth',2.0,...
      'defaultlinelinewidth',2.0,'defaultlinemarkersize',10.0)

figure(1)
plot(y1,y2, 'o');
hold on;
quiver(y1,y2, f1, f2)
axis([xmin-0.5 xmax+0.5 ymin-0.5 ymax+0.5])
