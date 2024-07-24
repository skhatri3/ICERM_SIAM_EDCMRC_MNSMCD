%main_channel

%Setting up Poiseuille flow velocity in channel and computing the
%force at given points 

%Developed by Shilpa Khatri and Ricardo Cortez  
%July 2024 

%Edited by Kristin Kurianskia and Brittany Leather

%clear all
%close all 

%% Parameters to set 

%setting the viscosity
mu = 1; 

%number of points on boundary where force is applied 
N = 200;

%choose blob
blob_num = 2;

%% Setting forces and computing velocity 

%channel
Lx = 10;
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

%top wall coordinates (x,y) = (y1_top,y2_top)
y1_top = s;
y2_top = ymax*ones(size(y1_top));

%bottom wall coordinates (x,y) = (y1_bot,y2_bot)
y1_bot = s;
y2_bot = ymin*ones(size(y1_bot));

%left-hand wall coordinates (x,y) = (y1_side,y2_side)
y2_side = (ymin+ds:ds:ymax-ds)';
y1_side = zeros(size(y2_side));

%coordiantes of boundary points
y1 = [y1_top; y1_bot; y1_side];
y2 = [y2_top; y2_bot; y2_side];

%velocity on boundary
%top
u1_top=zeros(size(y1_top));
u2_top=zeros(size(y1_top));
%bottom
u1_bot=zeros(size(y1_bot));
u2_bot=zeros(size(y1_bot));
%side
a=10;
u1_side=-a*(y2_side-ymin).*(y2_side-ymax);
u2_side=zeros(size(y1_side));

u1 = [u1_top; u1_bot; u1_side];
u2 = [u2_top; u2_bot; u2_side];

%computing the force 
f = RegStokeslets2D_velocitytoforce([y1,y2],[y1,y2],[u1,u2],ep,mu,blob_num);

%f is a force and to compare to exact solution need force density - divide
%by radius*dt 
f1 = f(:,1);
f2 = f(:,2);

%Points where velocity will be computed
x2 = y2_side;
x1 = 5*ones(size(x2));

%Compute velocity
u = RegStokeslets2D_forcetovelocity([y1, y2],f,[x1 x2],ep,mu,blob_num);


%% Compute velocity everywhere
%domain on which velocity is computed 
x1min = 0; 
x1max = xmax;
x2min = 0; 
x2max = ymax; 

%resolution for velocity grid
Nx1 = 80;  
Nx2 = 80;

%points on which velocity will be computed 
xx1 = linspace(x1min,x1max,Nx1);
xx2 = linspace(x2min,x2max,Nx2); 
[x1m,x2m] = ndgrid(xx1,xx2); 
x1 = x1m(:);
x2 = x2m(:);

%computing velocity 
uu = RegStokeslets2D_forcetovelocity([y1,y2],[f1,f2],[x1,x2],ep,mu,blob_num);
uu1 = uu(:,1);
uu2 = uu(:,2); 
u1m = reshape(uu1,size(xx1,2),size(xx2,2)); 
u2m = reshape(uu2,size(xx1,2),size(xx2,2));


%% Compute velocity at top boundary
%domain on which velocity is computed 
x1min = 0; 
x1max = xmax;
x2min = ymax; 
x2max = ymax; 

%resolution for velocity grid
Nx1 = 2000;  

%points on which velocity will be computed 
xx1 = linspace(x1min,x1max,Nx1)';
xx2 = ymax*ones(size(xx1));

%computing velocity 
uu = RegStokeslets2D_forcetovelocity([y1,y2],[f1,f2],[xx1,xx2],ep,mu,blob_num);

figure
plot(xx1, uu(:,2), '.-k')
%plot(xx1, sqrt(uu(:,1).^2+uu(:,2).^2), '.k')
ylim([-0.001,0.001])


%% Compute error
error1 = abs(u(:,1)-u1_side);
error2 = abs(u(:,2)-u2_side);

max(max(error1))/max(u1_side)
max(max(error2))

%figure; plot(x1, x2); hold on; quiver(u(:,1), u(:,2))

%% Plotting
% Plot channel with forces as quiver and velocity profile as dashed line
figure
plot(y1,y2,'k.')
hold on 
quiver(y1,y2,f1,f2,'r')

par = -a*(y2_side-ymin).*(y2_side-ymax);
plot(par, y2_side, 'b--')

axis equal
xlim([xmin,xmax])
ylim([ymin,ymax])

title('Forces and Computed Velocity')
xlabel('$x$', Interpreter='latex')
ylabel('$y$', Interpreter='latex')
ax = gca;
ax.FontSize = 20;


% Plot channel with velocity everywhere in channel
figure
plot(y1,y2,'k.')
surf(x1m, x2m, sqrt(u1m.^2+u2m.^2), 'EdgeColor', 'none')
view(2)

box on
colorbar 
axis equal 
grid off
xlim([xmin,xmax])
ylim([ymin,ymax])
title('Velocity Magnitude')
xlabel('$x$', Interpreter='latex')
ylabel('$y$', Interpreter='latex')