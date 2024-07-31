%main_channel

%Poiseuille flow velocity in a 2D channel

%Developed by Shilpa Khatri and Ricardo Cortez  
%July 2024 

%Edited by Kristin Kurianskia and Brittany Leather

clear all
close all 

%% Parameters to set 

%setting the viscosity
mu = 1; 

%number of points on boundary where force is applied 
N = 100;

%choose blob
blob_num = 1;

%resolution for velocity
Nx1 = 20*5;
Nx2 = 20;

% Setting forces and computing velocity 

%channel
Lx = 5;
Ly = 1;
xmin = -Lx;
xmax = Lx;
ymin = -Ly;
ymax = Ly;

%discretization step
ds = (xmax-xmin)/N;

%regularization parameter
ep = 1.5*ds;

%discretization of top and bottom
ds = (xmax-xmin)/N;
%stb = xmin+ds/2:ds:xmax-ds/2;
stb = xmin:ds:xmax;
stb = stb';

%discretization of left wall
%slr = ymin+ds/2:ds:ymax-ds/2;
slr = ymin+ds:ds:ymax-ds;
slr = slr';

%top wall coordinates (y1_top,y2_top)
y1_top = stb;
y2_top = ymax*ones(size(stb));

%bottom wall coordinates (y1_bot,y2_bot)
y1_bot = stb;
y2_bot = ymin*ones(size(stb));

%left wall coordinates (y1_left,y2_left)
y1_left = xmin*ones(size(slr));
y2_left = slr;

%right wall coordinates (y1_right,y2_right)
y1_right = xmax*ones(size(slr));
y2_right = slr;

%coordiantes of boundary points
y1 = [y1_top; y1_bot; y1_left; y1_right];
y2 = [y2_top; y2_bot; y2_left; y2_right];

%velocity on boundary
%top velocity
u1_top = zeros(size(stb));
u2_top = zeros(size(stb));

%bottom velocity
u1_bot = zeros(size(stb));
u2_bot = zeros(size(stb));

%left velocity
pL=1;
u_pois = -pL*(y2_left-ymin).*(y2_left-ymax);%/2/mu;

u1_left = u_pois;
u2_left = zeros(size(slr));

u1_right = u_pois;
u2_right = zeros(size(slr));

u1 = [u1_top; u1_bot; u1_left; u1_right];
u2 = [u2_top; u2_bot; u2_left; u2_right];

%computing the force 
f = RegStokeslets2D_velocitytoforce([y1,y2],[y1,y2],[u1,u2],ep,mu,blob_num);

%f is a force and to compare to exact solution need force density - divide
%by radius*dt 
f1 = f(:,1);
f2 = f(:,2);

%Points where velocity will be computed
xx1 = linspace(xmin,xmax,Nx1);
xx2 = linspace(ymin,ymax,Nx2);
[x1m, x2m] = ndgrid(xx1, xx2);
x1 = x1m(:);
x2 = x2m(:);

%computing velocity everywhere
ug = RegStokeslets2D_forcetovelocity([y1,y2],[f1,f2],[x1,x2],ep,mu,blob_num);
ug1 = ug(:,1);
ug2 = ug(:,2);
u1m = reshape(ug1,size(xx1,2),size(xx2,2));
u2m = reshape(ug2,size(xx1,2),size(xx2,2));
umag = sqrt(u1m.^2 + u2m.^2);

%compute pressure
p = RegStokeslets2D_forcetopressure([y1, y2],f,[x1 x2],ep,blob_num);
xpresh = xx1;
pmesh = reshape(p,size(xx1,2),size(xx2,2));
pmesh = pmesh-min(p);

%% Plotting
skip = 4;

figure(4)
plot(y1,y2,'k.')
hold on
quiver(y1,y2,u1,u2,'r')
axis equal

figure(5)
plot(y1,y2,'k.')
hold on
quiver(y1,y2,f1,f2,'r')
axis equal

figure(1)
plot(y1,y2,'k.')
hold on
pcolor(x1m,x2m,umag)
shading interp
quiver(y1,y2,f1,f2,'r')
axis equal
quiver(x1m(1:skip:end,1:skip:end),x2m(1:skip:end,1:skip:end),u1m(1:skip:end,1:skip:end),u2m(1:skip:end,1:skip:end),'k')
xlim([xmin,xmax])
ylim([ymin,ymax])
colorbar
title('Forces and Computed Velocity')