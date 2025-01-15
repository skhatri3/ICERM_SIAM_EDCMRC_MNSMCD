clear all 
close all

%% Parameters to set 

%setting the viscosity
mu = 1; 

%setting permeability 
Da = 0.3; 

%number of points on top (or bottom) wall 
N = 200;    

%resolution for velocity grid
Nx1 = 20*5;   
Nx2 = 20;

%blob choice - if you modify need to change beta below 
% phi = (2*d^4)/(pi*(r^2+d^2)^3) from 2021 paper
blob = 1;
% psi = 2*ep^4*(r^4 - 10*ep^2*r^2 + 5*ep^4)/(Pi*(r^2 + ep^2)^5) from 2021 paper 
%blob = 2; 

%channel
Lx = 10;
Ly = 1;
xmin = 0;
xmax = Lx;
ymin = -Ly;
ymax = Ly;

%discretization of channel
ds = (xmax-xmin)/N;

%regularization parameter
ep = 0.3; %1.5*ds; %0.15; %3*ds; %0.0224;%1*ds;

%discretization of top and bottom:
%stb = xmin+ds/2:ds:xmax-ds/2;
stb = xmin:ds:xmax; 
stb = stb';

st = 0:ds:10; 
st = st'; 

%discretization of left and right 
%slr = ymin+ds/2:ds:ymax-ds/2;
slr = ymin+ds:ds:ymax-ds;
slr = slr'; 

%top wall coordinates (y1_top,y2_top)
y1_top = st;
y2_top = ymax*ones(size(st));

%unit normals for top: 
normals_top = zeros(length(st),2);
normals_top(:,2) = 1;

%bottom wall coordinates (y1_bot,y2_bot)
y1_bot = stb;
y2_bot = ymin*ones(size(stb));

%unit normals for bottom
normals_bot = zeros(length(stb),2); 
normals_bot(:,2) = -1;

%left wall coordinates (y1_left,y2_left)
y1_left = xmin*ones(size(slr));
y2_left = slr;

%unit normals for left  wall 
normals_left = zeros(length(slr),2);
normals_left(:,1) = -1;

%right wall coordinates (y1_right,y2_right)
y1_right = xmax*ones(size(slr));
y2_right = slr;

%unit normals for right wall 
normals_right = zeros(length(slr),2);
normals_right(:,1) = 1;

%combining boundary points 
y1 = [y1_top; y1_bot]; %; y1_left; y1_right];
y2 = [y2_top; y2_bot]; %; y2_left; y2_right];

%y1 = [y1_bot; y1_left; y1_right];
%y2 = [y2_bot; y2_left; y2_right];

%y1 = [y1_top; y1_bot; y1_left];
%y2 = [y2_top; y2_bot; y2_left];

%permeable region in y1 and y2 
I = [1:length([ y1_top])]';
    
%combining normals 
normals=[normals_top; normals_bot]; %; normals_left; normals_right];

%normals=[normals_top; normals_bot; normals_left];

%beta vector for function;
%beta = -0.*sign([y1_top; y1_bot]).*4*ep*alpha.*ones(length([ y1_top; y1_bot ]),1)/3;  
%beta = Da.*ones(length([ y1_top; y1_bot ]),1);
%beta = 10*ones(length([ y1_top; y1_bot ]),1); 
%beta = Da.*ones(length([ y1_top]),1)/500;
%beta = 0.00045*ones(length([ y1_top]),1)./ep; 
beta = 0.003*ones(length([ y1_top]),1); 
%beta = [beta; zeros(length([ y1_bot; y1_left; y1_right ]),1) ];
beta = [beta; zeros(length([ y1_bot]),1)];
%beta = [beta; zeros(length(y1_left),1) ];
%beta = beta./ds; 

%velocity on boundary

%top
u1_top = zeros(size(st));
u2_top = zeros(size(st));

%bottom
u1_bot = zeros(size(stb));
u2_bot = ones(size(stb));

%left: poiseulle flow
u1_left = zeros(size(slr)); 
u2_left = ones(size(slr));

%right: poiseulle flow
u1_right = zeros(size(slr)); 
u2_right = ones(size(slr));

%combining velocities 
u1 = [u1_top; u1_bot]; %; u1_left; u1_right];
u2 = [u2_top; u2_bot]; %; u2_left; u2_right];

%u1 = [u1_bot; u1_left; u1_right];
%u2 = [u2_bot; u2_left; u2_right];

%u1 = [u1_top; u1_bot; u1_left];
%u2 = [u2_top; u2_bot; u2_left];

%compute g by setting boundary conditions (Eq 19 in 2021 paper) 
g = RegStokeslets2D_velocityto_gforce_permeable ([y1,y2],[y1,y2],...
    [u1,u2], ep, mu, blob, I, beta, normals);

%find velocity in permeable region:
[u_perm]=RegStokeslets2D_permeable_gtovelocity([y1,y2], g, [y1(I),y2(I)],...
    ep, mu, blob, beta, normals);
u_perm1 = u_perm(:,1);
u_perm2 = u_perm(:,2); 

%put in the computed velocities for the permeable part
%first save initial velocities 
u1_init = u1; 
u2_init = u2; 
%u1 = u1*0; 
%u2 = u2*0;
u1(I) = u_perm1; 
u2(I) = u_perm2; 

%computing the force with boundary velocities only using Stokeslets 
f = RegStokeslets2D_velocitytoforce([y1,y2],[y1,y2],[u1,u2],ep,mu,blob);
f1 = f(:,1);
f2 = f(:,2);  

%points on which velocity will be computed 
xx1 = linspace(xmin,xmax,Nx1);
xx2 = linspace(ymin,ymax,Nx2); 
[x1m,x2m] = ndgrid(xx1,xx2); 
x1 = x1m(:);
x2 = x2m(:);

%computing velocity 
ug = RegStokeslets2D_forcetovelocity([y1,y2],[f1,f2],[x1,x2],ep,mu,blob);
ug1 = ug(:,1);
ug2 = ug(:,2); 
u1m = reshape(ug1,size(xx1,2),size(xx2,2)); 
u2m = reshape(ug2,size(xx1,2),size(xx2,2));
umag = sqrt(u1m.^2 + u2m.^2);

% %exact solution 
% [uexact1, uexact2] = channel_permeable_infinite_exact_solution(x1,x2);
% uexact1m = reshape(uexact1,size(xx1,2),size(xx2,2)); 
% uexact2m = reshape(uexact2,size(xx1,2),size(xx2,2)); 
% 
% %computing error 
% error1 = abs(u1m-uexact1m);
% error2 = abs(u2m-uexact2m); 
% errormag = sqrt(error1.^2 + error2.^2); 
% 
% %prints the max error 
% fprintf('maximum error in u1: %d \n',max(max(error1)));
% fprintf('maximum error in u2: %d \n',max(max(error2)));

%% Plotting figures 
skip = 2; %for quiver plots 
set(0,'defaultaxesfontsize',20,'defaultaxeslinewidth',2.0,...
      'defaultlinelinewidth',2.0,'defaultlinemarkersize',10.0)

figure(6)
plot(y1,y2,'k.')
hold on
quiver(y1,y2,normals(:,1),normals(:,2),'b')
axis equal 

figure(2)
plot(y1,y2,'k.')
hold on
quiver(y1,y2,u1_init,u2_init,'r')
axis equal 

figure(3) 
plot(y1,y2,'k.')
hold on
quiver(y1,y2,g(:,1),g(:,2),'r')
axis equal

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

% figure(6)
% plot(y1,y2,'k.')
% hold on 
% pcolor(x1m,x2m,log10(errormag)+eps)
% shading interp
% colorbar

