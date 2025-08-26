% Flow through permeable channel using exact solution on boundaries and
% left/right 

%Adapted from main_channel by Kristin Kurianski Aug 2025
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
[u1_top,u2_top,~] = permeablechannelexact(y1_top,y2_top);

%bottom velocity
[u1_bot,u2_bot,~] = permeablechannelexact(y1_bot,y2_bot);

%left velocity
% pL=1;
% u_pois = -pL*(y2_left-ymin).*(y2_left-ymax);%/2/mu;
%u1_left = u_pois;
%u2_left = zeros(size(slr));
%u1_right = u_pois;
%u2_right = zeros(size(slr));
[u1_left,u2_left,~] = permeablechannelexact(y1_left,y2_left);
[u1_right,u2_right,~] = permeablechannelexact(y1_right,y2_right);

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

%exact pressure
[~,~,p_exact] = permeablechannelexact(x1m,x2m);

%% Plotting
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
figure(5)
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