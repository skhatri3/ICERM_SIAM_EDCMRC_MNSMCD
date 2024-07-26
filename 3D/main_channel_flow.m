% Setting velocities on a cylindrical channel of radius a and computing the
% forces at given points and then the velocity interior to the channel 
% Velocity on channels is set to 0 on side walls, gamma*(a^2-r^2) on the
% inlet and on the outlet 

% Developed for ICERM EDCMRC MNSMCD 
% July 2024 

function [l2error1,l2error2,l2error3,maxerror1,maxerror2,maxerror3]= main_channel_flow(Nt,c_ep)
%inputs 
%Nt: approx number of points to discretize cylinder in the radial direction 
%Number of points in the length is then automatically determined 
%c_ep: constant multiplying discretization for blob size 

%setting the viscosity
mu = 1; 

%radius of cylinder 
a = 1; 

%length of cylinder 
L = 5; 

% 2d surface (x1-x3 plane) on which velocity is computed 
x1min = 0; 
x1max = L; 
x2fixed = 0; 
x3min = -a; 
x3max = a; 

%resolution for velocity grid
Nx1 = 200;  
Nx3 = 200; 

%discretization of cylinder 
[y1,y2,y3,darea,Ntube,Ncap] = cylinder_surface_twocaps(Nt,L,a);

%regularization parameter
%c_ep = 1.5; 
ep = c_ep*mean(sqrt(darea))^(0.5); 

%total  number of points on surface 
Npts = Ntube + 2*Ncap; 

%velocity on tube part of cylinder is set to 0  
us1 = 0*ones(Ntube,1); 
us2 = 0*ones(Npts,1); 
us3 = 0*ones(Npts,1);

for i = Ntube+1:Ntube+2*Ncap 
    rcap2 = y2(i)^2 + y3(i)^2; 
    us1(i) = 1 - rcap2/a/a; 
end 

%computing force on surface of sphere 
fs = RegStokeslets3D_velocitytoforce([y1,y2,y3],[y1,y2,y3],[us1,us2,us3],ep,mu);
f1 = fs(:,1);
f2 = fs(:,2); 
f3 = fs(:,3); 

%points on which velocity will be computed 
xx1 = linspace(x1min,x1max,Nx1);
xx3 = linspace(x3min,x3max,Nx3); 
[x1m,x3m] = ndgrid(xx1,xx3); 
x2m = x1m*0 + x2fixed; 
x1 = x1m(:);
x2 = x2m(:);
x3 = x3m(:);

%discretization of plane on which velocity will be computed 
h1 = (x1max-x1min)/Nx1; 
h3 = (x3max-x3min)/Nx3; 

%area of plane on which velocity will be computed 
pa = (x3max-x3min)*(x1max-x1min);

%computing velocity on 2D surface
u = RegStokeslets3D_forcetovelocity([y1,y2,y3],[f1,f2,f3],[x1,x2,x3],ep,mu);
u1 = u(:,1);
u2 = u(:,2); 
u3 = u(:,3); 
u1m = reshape(u1,size(xx1,2),size(xx3,2)); 
u2m = reshape(u2,size(xx1,2),size(xx3,2));
u3m = reshape(u3,size(xx1,2),size(xx3,2));

%% Computing error on x1-x3 plane

%exact solution on x1-x3 plane - assumes only interior to cylinder 

for i = 1:length(xx1)

    for j = 1:length(xx3)

      r2 = x2m(i,j).^2 + x3m(i,j).^2; %radius 
    
      uexact1(i,j) = (1-r2/a/a); 
      uexact2(i,j) = 0; 
      uexact3(i,j) = 0;

    end
end

%computing error 
error1 = abs(u1m-uexact1);
error2 = abs(u2m-uexact2); 
error3 = abs(u3m-uexact3);
errormag = sqrt(error1.^2 + error2.^2 + error3.^2); 
l2error1 = sqrt(sum(sum(error1.^2*h1*h3))/pa); 
l2error2 = sqrt(sum(sum(error2.^2*h1*h3))/pa); 
l2error3 = sqrt(sum(sum(error3.^2*h1*h3))/pa); 
maxerror1 = max(max(error1)); 
maxerror2 = max(max(error2)); 
maxerror3 = max(max(error3));

fprintf('Error in velocity on 2D plane \n')

%prints the max error 
fprintf('maximum error in u1: %d \n',maxerror1);
fprintf('maximum error in u2: %d \n',maxerror2);
fprintf('maximum error in u3: %d \n',maxerror3);

%prints the l2 error 
fprintf('l2 error in u1: %d \n',l2error1);
fprintf('l2 error in u2: %d \n',l2error2);
fprintf('l2 error in u3: %d \n',l2error3);


%% Plotting figures 
skip = 8; %for quiver plots 
set(0,'defaultaxesfontsize',20,'defaultaxeslinewidth',2.0,...
      'defaultlinelinewidth',2.0,'defaultlinemarkersize',10.0)


%plotting on 2D surface in domain 

figure(1)
plot3(y1,y2,y3,'.')
hold on 
quiver3(y1,y2,y3,f1,f2,f3)
axis equal 
%splot = surf(x1m,x2m,x3m,u1m);
%splot.EdgeColor = 'none';
%colorbar
%quiver3(x1m(1:skip:end,1:skip:end),x2m(1:skip:end,1:skip:end),x3m(1:skip:end,1:skip:end),u1m(1:skip:end,1:skip:end),u2m(1:skip:end,1:skip:end),u3m(1:skip:end,1:skip:end),'k')
%view(0,0)
%title('Numerical Solution')

figure(2) 
plot3(y1,y2,y3,'.')
hold on 
%quiver3(y1,y2,y3,f1,f2,f3)
axis equal 
quiver3(x1m(1:skip:end,1:skip:end),x2m(1:skip:end,1:skip:end),x3m(1:skip:end,1:skip:end),u1m(1:skip:end,1:skip:end),u2m(1:skip:end,1:skip:end),u3m(1:skip:end,1:skip:end),'k')

figure(3)
pcolor(x1m,x3m,u1m)
shading interp 
hold on 
quiver(x1m(1:skip:end,1:skip:end),x3m(1:skip:end,1:skip:end),u1m(1:skip:end,1:skip:end),u3m(1:skip:end,1:skip:end),'k')
axis equal
colorbar 

figure(4)
pcolor(x1m,x3m,u2m)
shading interp 
hold on 
quiver(x1m(1:skip:end,1:skip:end),x3m(1:skip:end,1:skip:end),u1m(1:skip:end,1:skip:end),u3m(1:skip:end,1:skip:end),'k')
axis equal
colorbar 

figure(5)
pcolor(x1m,x3m,u3m)
shading interp 
hold on 
quiver(x1m(1:skip:end,1:skip:end),x3m(1:skip:end,1:skip:end),u1m(1:skip:end,1:skip:end),u3m(1:skip:end,1:skip:end),'k')
axis equal
colorbar 

figure(6)
pcolor(x1m,x3m,log10(errormag)+eps)
shading interp
colorbar
axis equal 
clim([-5,0])
title('Error')

