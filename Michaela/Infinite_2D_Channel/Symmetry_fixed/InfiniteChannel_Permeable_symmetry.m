% Flow through permeable channel using exact solutions on left and right
% boundaries from Bernardi et al 2023.

% Adapted from main_channel by Kristin Kurianski and Michaela Kubacki Fall
% 2025. Original code: main_channel for Poiseuille flow velocity in a 2D
% channel was developed by Shilpa Khatri and Ricardo Cortez July 2024

% Recent Updates:  
%   - 12/2: updated all function files to consistently rescale the 
%     Stokeslet and source double matrices
%   - 12/3: updated this file and all function files to incorporate quad.
%     weight vector, wt, so that we have symmetry.  Originally we were
%     using same discretization width on sides and top/bottom, which meant
%     that the source points on sides were not symmetrically spaced.
%   - Just guessing regarding the relationship between beta with Da and ep.
%     Currently using beta = C*Da*ep and manually guessing C. This was
%     similar to the relationship found in previous applications.  However,
%     note that if you change Da the relationship doesn't hold!


clear all
close all

%% Parameters to set

% Model Parameters
Da = 0.4; % Darcy Number
mu = 1; % Viscosity
blob_num = 2; % blob choice

% Number of source and target points
N = 50; % Number of source points (along top  and bottom boundaries)
Nx1 = 40*5; % Number of target points in x direction for full channel calc
Nx2 = 40; % Number of target points in y direction for full channel calc

% Setting up channel

%                 1 - top (perm)
%         +""""""""""""""""""""""""""""""+
%         |                              |
%  3-left |                              | 4-right
% (known) |                              | (known)
%         |                              |                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     |   
%         +""""""""""""""""""""""""""""""+
%                2 - bottom (perm)

% Channel Geometry
Lx = pi;
Ly = 1;
xmin = -Lx;
xmax = Lx;
ymin = -Ly;
ymax = Ly;

% Discretization step
%ds = (xmax-xmin)/N;
ds_y = (ymax - ymin)/N;
ds_x = (xmax - xmin)/( ceil((xmax - xmin)/ds_y));
% Regularization parameter
ep = ds_y;

% Discretize Channel Walls (using midpoint to avoid the corners)
% top and bottom
stb = xmin+ds_x/2:ds_x:xmax-ds_x/2;
stb = stb';
% left and right
slr = ymin+ds_y/2:ds_y:ymax-ds_y/2;
slr = slr';

% Define coordinates on each wall (x-coord,y-coord)
y1_top = stb;   y2_top = ymax*ones(size(stb)); % top wall coordinates (y1_top,y2_top)
y1_bot = stb;   y2_bot = ymin*ones(size(stb)); % bottom wall coordinates (y1_bot,y2_bot)
y1_left = xmin*ones(size(slr));   y2_left = slr; % left wall coordinates (y1_left,y2_left)
y1_right = xmax*ones(size(slr));   y2_right = slr; % right wall coordinates (y1_right,y2_right)
y1 = [y1_top; y1_bot; y1_left; y1_right]; %x-coordinates of all boundary points
y2 = [y2_top; y2_bot; y2_left; y2_right]; %y-coordinates of all boundary points

% Define weights cooresponding to wall coordinates
wt = [ds_x*ones(size(y1_top)); ds_x*ones(size(y1_bot)); ds_y*ones(size(y1_left)); ds_y*ones(size(y1_right))];

y1_tb = [y1_top; y1_bot]; % Top and bottom x coordinates
y2_tb = [y2_top; y2_bot]; % Top and bottom y coordinates

% Exact velocities on boundaries using true solution (Bernardi, 2023)
[u1_top_exact,u2_top_exact,p_top_exact] = permeablechannelexact(y1_top,y2_top,Da); %top velocity
[u1_bot_exact,u2_bot_exact,p_bot_exact] = permeablechannelexact(y1_bot,y2_bot,Da); %bottom velocity
[u1_left_exact,u2_left_exact,p_left_exact] = permeablechannelexact(y1_left,y2_left,Da); % left velocity
[u1_right_exact,u2_right_exact,p_right_exact] = permeablechannelexact(y1_right,y2_right,Da); %right velocity

u1_bd_exact = [u1_top_exact; u1_bot_exact; u1_left_exact; u1_right_exact]; %x-coordinates of all boundary velocities
u2_bd_exact  = [u2_top_exact; u2_bot_exact; u2_left_exact; u2_right_exact]; %y-coordinates of all boundary velocities
p_bd_exact = [p_top_exact; p_bot_exact; p_left_exact; p_right_exact]; % pressure on boundary

u1_tb_exact = [u1_top_exact; u1_bot_exact];
u2_tb_exact = [u2_top_exact; u2_bot_exact];
p_tb_exact = [p_top_exact; p_bot_exact];

% Quiver plot of Exact velocities on boundaries
figure
plot(y1,y2,'k.')
hold on
quiver(y1,y2,u1_bd_exact,u2_bd_exact,'r')
title('Exact Velocities on Boundaries')
set(gca,'Fontsize',20)
axis equal
axis([-4,4,-1.5,1.5])

% Define normal vectors
normals_top = zeros(length(stb),2); % unit normals for top:
normals_top(:,2) = 1;
normals_bot = zeros(length(stb),2); % unit normals for bottom
normals_bot(:,2) = -1;
normals_left = zeros(length(slr),2); % unit normals for left  wall 
normals_left(:,1) = -1;
normals_right = zeros(length(slr),2); % unit normals for right wall 
normals_right(:,1) = 1;
% normals on full boundary
normals = [normals_top; normals_bot; normals_left; normals_right];
% normals for top and bottom only
normals_tb = [normals_top; normals_bot];

% indices of permeable region (top and bottom boundaries)
idx = (1:2*length( y1_top ))';

%% Stokeslets Only
%
% Use Stokeslets and exact velocities on top and bottom to find f on Gamma
% Given the exact velocities on boundary of channel (from the true
% solution), find the forces on boundary of channel, then use the boundary
% forces to find the velocities throughout the channel

% All four sides of channel
ftemp = RegStokeslets2D_velocitytoforce([y1,y2], [y1,y2], [u1_bd_exact,u2_bd_exact], ep, mu, blob_num,wt);

% Quiver Plot of Forces
% figure
% plot(y1,y2,'k.')
% hold on
% quiver(y1,y2,ftemp(:,1),ftemp(:,2),'b')
% title('Forces on Boundaries via Stokeslets Only')
% set(gca,'Fontsize',20)
% axis equal
% axis([-4,4,-1.5,1.5])

%Points where velocity will be computed everywhere else
xx1 = linspace(xmin+ds_x/2,xmax-ds_x/2,Nx1);
xx2 = linspace(ymin+ds_y/2,ymax-ds_y/2,Nx2);
[x1m, x2m] = ndgrid(xx1, xx2);
x1 = x1m(:); %x-coords of points where computing velocity
x2 = x2m(:); %y-coords of points where computing velocity

% From boundary forces, calculate channel velocities
ug = RegStokeslets2D_forcetovelocity([y1,y2],ftemp,[x1,x2],ep,mu,blob_num,wt);

ug1 = ug(:,1);
ug2 = ug(:,2);
u1m = reshape(ug1,size(xx1,2),size(xx2,2)); %x-coords of computed velocities
u2m = reshape(ug2,size(xx1,2),size(xx2,2)); %y-coords of computed velocities

% Plot streamlines comparing exact vs computed velocities using Stokeslets
figure
plot(y1, y2, 'k.')
hold on
% get exact solution everywhere
[x1gg, x2gg] = meshgrid(xx1, xx2);
[uexact,vexact,~] = permeablechannelexact(x1gg,x2gg,Da);

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

title('Computed (blue) and exact (yellow) streamlines - Stokeslets Only')
set(gca, 'FontSize', 20)
xlabel('$x$', 'Interpreter', 'latex')
ylabel('$y$', 'Interpreter', 'latex')
axis equal
axis([-4,4,-1.5,1.5])


%% Permeable Channel Problem
% Now we utilize the Stokeslets + Source Doublets approach to simulate the
% permeable channel flow

% Approach: Use Darcy number times epsilon for beta
beta = zeros (length(y1),1);
beta(idx) = Da*ep*0.039; % right now guessing beta = Da*ep*constant

% 4-Step Permeable Channel Process: Input velocities are the known (exact)
% velocities on left and right, and zero for unknown velocities on top and
% bottom (where it's permeable).

% Velocity placeholders for top and bottom boundaries
u1_top = zeros(size(y1_top)); %top velocity (no contribution from stokeslets)
u2_top = zeros(size(y1_top)); %top velocity (no contribution from stokeslets)
u1_bot = zeros(size(y1_top)); %top velocity (no contribution from stokeslets)
u2_bot = zeros(size(y1_top)); %top velocity (no contribution from stokeslets)

u1_tb = [u1_top; u1_bot]; % top and bottom velocities, x-coordinates
u2_tb = [u2_top; u2_bot]; % top and bottom velocities, y-coordinates

u1 = [u1_top; u1_bot; u1_left_exact; u1_right_exact]; %x-coordinates of all boundary velocities
u2 = [u2_top; u2_bot; u2_left_exact; u2_right_exact]; %y-coordinates of all boundary velocities

% Step 1: Find g-force distribution due to St+SD on solid walls and SD (no
% Stokeslets) on permeable walls

[g] = RegStokeslets2D_velocityto_gforce_permeable([y1,y2],[y1,y2],...
    [u1,u2], ep, mu, blob_num, idx, beta, normals, wt);

% Step 2: Use g on full boundary to find missing velocities on permeable
% walls

y1b=y1(idx);
y2b=y2(idx);

[u_perm] = RegStokeslets2D_gtovelocity([y1,y2], g, [y1_tb,y2_tb],...
    ep, mu, blob_num, beta, normals, wt);
u_perm1 = u_perm(:,1);
u_perm2 = u_perm(:,2);

% Replace the placeholder zeros in u1,u2 with the velocities
u1(idx) = u_perm1;
u2(idx) = u_perm2;

% Quiver plot of recomputed velocities on boundaries
figure
subplot(2,1,1);
hold on;
plot(y1,y2,'k.')
hold on
quiver(y1,y2,u1,u2,'b')
%quiver([y1_left; y1_right],[y2_left; y2_right],u1(203:end),u2(203:end),'b')
title('Computed Velocities on Boundaries after using g')
set(gca,'Fontsize',14)
hold off;
axis equal;
axis([-4,4,-1.5,1.5])
subplot(2,1,2);
plot(y1,y2,'k.')
hold on
quiver(y1,y2,u1_bd_exact, u2_bd_exact,'r')
hold off
title('Exact Boundary Velocities')
set(gca,'Fontsize',14)
axis equal
axis([-4,4,-1.5,1.5])


% Now that the missing velocities are recovered, steps 3-4 are the regular
% "Stokeslet" problem, given boundary velocities, find boundary forces,
% then find channel velocities.

% Step 3: Use boundary velocities to solve for forces on boundary
f = RegStokeslets2D_velocitytoforce([y1,y2], [y1,y2], [u1,u2], ep, mu, blob_num, wt);
f1 = f(:,1);
f2 = f(:,2);

% Step 4: compute velocity everywhere using the computed forces

ug = RegStokeslets2D_forcetovelocity([y1,y2],[f1,f2],[x1,x2],ep,mu,blob_num, wt);
ug1 = ug(:,1);
ug2 = ug(:,2);
u1m = reshape(ug1,size(xx1,2),size(xx2,2)); %x-coords of computed velocities
u2m = reshape(ug2,size(xx1,2),size(xx2,2)); %y-coords of computed velocities
umag = sqrt(u1m.^2 + u2m.^2);

% Exact solutions
[u1_exact,u2_exact,p_exact] = permeablechannelexact(x1m,x2m,Da);
umag_exact = sqrt(u1_exact.^2 + u2_exact.^2);

% Plot streamlines
figure
plot(y1, y2, 'k.')
hold on
% get exact solution everywhere
[x1gg, x2gg] = meshgrid(xx1, xx2);
[uexact,vexact,~] = permeablechannelexact(x1gg,x2gg,Da);

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

title('Computed (blue) and exact (yellow) streamlines - Permeable Channel')
set(gca, 'FontSize', 20)
xlabel('$x$', 'Interpreter', 'latex')
ylabel('$y$', 'Interpreter', 'latex')
axis([-4,4,-1.5,1.5])

%% Exact solutions
function [ug,vg,pg] = permeablechannelexact(x,y,Da)

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