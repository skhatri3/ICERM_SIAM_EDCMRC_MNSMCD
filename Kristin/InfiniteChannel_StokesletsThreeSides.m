% Flow through permeable channel using exact solution on three sides of
% channel with doublets on the top

%Adapted from main_channel by Kristin Kurianski Sep 2025
%main_channel for Poiseuille flow velocity in a 2D channel
%Developed by Shilpa Khatri and Ricardo Cortez
%July 2024

clear all
close all

%% Parameters to set
%setting the viscosity
mu = 1;

%number of points on boundary where force is applied
N = 100;

%choose blob
blob_num = 2;

%resolution for velocity
% Computing error using velocity in x-direction
Nx1 = 20*5;
Nx2 = 20;

% Setting forces and computing velocity

%channel
Lx = pi;
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
%stb = xmin+ds/2:ds:xmax-ds/2;
stb = xmin:ds:xmax;
stb = stb';
%discretization of left and right
%slr = ymin+ds/2:ds:ymax-ds/2;
slr = ymin+ds:ds:ymax-ds;
slr = slr';

% define points on each wall (x-coord,y-coord)
y1_top = stb;   y2_top = ymax*ones(size(stb)); %top wall coordinates (y1_top,y2_top)
y1_bot = stb;   y2_bot = ymin*ones(size(stb)); %bottom wall coordinates (y1_bot,y2_bot)
y1_left = xmin*ones(size(slr));   y2_left = slr; %left wall coordinates (y1_left,y2_left)
y1_right = xmax*ones(size(slr));   y2_right = slr; %right wall coordinates (y1_right,y2_right)
y1 = [y1_top; y1_bot; y1_left; y1_right]; %x-coordiantes of all boundary points
y2 = [y2_top; y2_bot; y2_left; y2_right]; %y-coordiantes of all boundary points

%set velocity on boundary using exact solution (Bernardi, 2023)
%[u1_top,u2_top,~] = permeablechannelexact(y1_top,y2_top); %top velocity
u1_top = zeros(size(y1_top)); %top velocity (no contribution from stokeslets)
u2_top = zeros(size(y1_top)); %top velocity (no contribution from stokeslets)
[u1_bot,u2_bot,~] = permeablechannelexact(y1_bot,y2_bot); %bottom velocity
[u1_left,u2_left,~] = permeablechannelexact(y1_left,y2_left); % left velocity
[u1_right,u2_right,~] = permeablechannelexact(y1_right,y2_right); %right velocity
u1 = [u1_top; u1_bot; u1_left; u1_right]; %x-coordiantes of all boundary velocities
u2 = [u2_top; u2_bot; u2_left; u2_right]; %y-coordinates of all boundary velocities

%define normal vectors for permeable region
%unit normals for top:
normals_top = zeros(length(stb),2);
normals_top(:,2) = 1;
%unit normals for bottom
normals_bot = zeros(length(stb),2); 
normals_bot(:,2) = -1;
%unit normals for left  wall 
normals_left = zeros(length(slr),2);
normals_left(:,1) = -1;
%unit normals for right wall 
normals_right = zeros(length(slr),2);
normals_right(:,1) = 1;
% all normals along boundaries
normals = [normals_top; normals_bot; normals_left; normals_right];

% indices of permeable region in y1 and y2
idx = (1:length( y1_top ))';

% set beta
beta = zeros(length(y1), 1);
beta(idx) = (1*ep);

%compute g by setting boundary conditions
g = RegStokeslets2D_velocityto_gforce_permeable ([y1,y2],[y1,y2],...
    [u1,u2], ep, mu, blob_num, idx, beta, normals);

%find velocity in permeable region:
y1b=y1(idx);
y2b=y2(idx);
[u_perm] = RegStokeslets2D_permeable_gtovelocity([y1,y2], g, [y1(idx),y2(idx)],...
    ep, mu, blob_num, beta, normals);
u_perm1 = u_perm(:,1);
u_perm2 = u_perm(:,2);

%put in the computed velocities for the permeable part
%first save initial velocities
u1_init = u1;
u2_init = u2;
u1(idx) = u_perm1;
u2(idx) = u_perm2;

%computing the force with boundary velocities only using Stokeslets
f = RegStokeslets2D_velocitytoforce([y1,y2], [y1,y2], [u1,u2], ep, mu, blob_num);
%f is a force and to compare to exact solution need force density - divide
%by radius*dt
f1 = f(:,1);
f2 = f(:,2);

%Points where velocity will be computed everywhere else
xx1 = linspace(xmin,xmax,Nx1);
xx2 = linspace(ymin,ymax,Nx2);
[x1m, x2m] = ndgrid(xx1, xx2);
x1 = x1m(:); %x-coords of points where computing velocity
x2 = x2m(:); %y-coords of points where computing velocity

%computing velocity everywhere using the computed forces
ug = RegStokeslets2D_forcetovelocity([y1,y2],[f1,f2],[x1,x2],ep,mu,blob_num);
ug1 = ug(:,1);
ug2 = ug(:,2);
u1m = reshape(ug1,size(xx1,2),size(xx2,2)); %x-coords of computed velocities
u2m = reshape(ug2,size(xx1,2),size(xx2,2)); %y-coords of computed velocities
umag = sqrt(u1m.^2 + u2m.^2);

%exact solutions
[u1_exact,u2_exact,p_exact] = permeablechannelexact(x1m,x2m);
umag_exact = sqrt(u1_exact.^2 + u2_exact.^2);

% compute magnitude of residual along center of velocity field (x=0)
% error_centermag = zeros(1,length(N));
% middle_ind = floor(Nx1/2);
% error_centermag = norm(umag(50,:) - umag_exact(50,:));




%% Plotting
% residual of velocity in center of x-interval
% figure(6)
% scatter(N, error_centermag, 50, 'filled', 'linewidth', 2);
% title('Magnitude of residual of velocity in center of x-interval')
% xlabel('Number of boundary points where forces applied')
% ylabel('$\|u(0,y) - u_{exact}(0,y)\|_2$', 'interpreter', 'latex')
% ax=gca;
% ax.FontSize=14;
% 
% figure(7)
% semilogy(N,error_centermag, 'bo', 'linewidth',2, 'markersize', 6, 'markerfacecolor', 'b')
% title('Log plot of agnitude of residual of velocity in center of x-interval')
% xlabel('Number of boundary points where forces applied')
% ylabel('$\log(\|u(0,y) - u_{exact}(0,y)\|_2)$', 'interpreter', 'latex')
% ax=gca;
% ax.FontSize=14;



% Velocity and force plots
skip = 1;

% Quiver plot of exact velocities on boundaries
figure(1)
plot(y1,y2,'k.')
hold on
quiver(y1,y2,u1,u2,'r')
axis equal

% Quiver plot of computed forces on boundaries
figure(2)
plot(y1,y2,'k.')
hold on
quiver(y1,y2,f1,f2,'r')
axis equal

% pcolor plot of velocities everywhere with quiver of computed forces on boundaries
figure(3)
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

% Plot exact pressure and exact boundary velocities
figure(4)
hold on
contourf(x1m,x2m,p_exact, 'EdgeColor', 'none')
colormap sky
colorbar
% boundary points
plot(y1(1:skip:end), y2(1:skip:end), 'k.', 'LineWidth', 1)
% quiver of exact velocities on boundaries
quiver(y1(1:skip:end), y2(1:skip:end), u1(1:skip:end, 1:skip:end), u2(1:skip:end, 1:skip:end), 1, 'k-', 'LineWidth', 1)
% plot attributes
xlabel('$x$', 'Interpreter', 'latex')
ylabel('$y$', 'Interpreter', 'latex')
axis equal,axis([-3.5 3.5 -1.5 1.5])
box on
legend('Exact pressure', '', 'Exact boundary velocity')
ax = gca;
ax.FontSize = 14;
hold off

% Plot streamlines
figure(6)
plot(y1, y2, 'k.')
hold on
% get exact solution everywhere
[x1gg, x2gg] = meshgrid(xx1, xx2);
[uexact,vexact,~] = permeablechannelexact(x1gg,x2gg);

% Computed streamlines
hh1_comp = streamline(x1gg, x2gg, u1m', u2m', x1gg(1:end, 1), x2gg(1:end, 1)); % streamlines starting at left wall
hh2_comp = streamline(x1gg, x2gg, -u1m', -u2m', x1gg(1:end, end), x2gg(1:end, end)); % streamlines starting at right wall
% Exact streamlines
hh1_exact = streamline(x1gg, x2gg, uexact, vexact, x1gg(1:end, 1), x2gg(1:end, 1)); % streamlines starting at left wall
hh2_exact = streamline(x1gg, x2gg, -uexact, -vexact, x1gg(1:end, end), x2gg(1:end, end)); % streamlines starting at right wall


strmLW = 2;
RGB = orderedcolors("gem");
set(hh1_comp, 'Color', RGB(1,:), 'linewidth', strmLW);
set(hh2_comp, 'Color', RGB(1,:), 'linewidth', strmLW);
set(hh1_exact, 'Color', RGB(3,:), 'LineStyle','--', 'linewidth', strmLW);
set(hh2_exact, 'Color', RGB(3,:), 'LineStyle','--', 'linewidth', strmLW);


title('Computed (blue) and exact (yellow) streamlines')
set(gca, 'FontSize', 20)
xlabel('$x$', 'Interpreter', 'latex')
ylabel('$y$', 'Interpreter', 'latex')
axis equal;

function [ug,vg,pg] = permeablechannelexact(x,y)
%set Darcy number
Da = 0.3;

% Newton's method to find lambda1
L1=sqrt(2*Da);
f=@(L) tan(L)^2-tan(L)/L+1-2*Da;
L1 = fzero(f,L1);
%disp(['lambda_1 = ',num2str(L1)])

% Compute g and g' from FB paper
g  = cot(L1)*(Da-0.5)   *sin(L1*y)+0.5*y.*cos(L1*y);
gp = cot(L1)*(Da-0.5)*L1*cos(L1*y)+0.5*cos(L1*y)-0.5*L1*y.*sin(L1*y);

% Compute exact solution using g and g'
C = -4;
pg = C*sinh(L1*x) .* cos(L1*y);
ug = -C/L1 * cosh(L1*x) .* gp;
vg =  C    * sinh(L1*x) .* g ;

end