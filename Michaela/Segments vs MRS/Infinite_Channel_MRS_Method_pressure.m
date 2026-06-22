% Flow through permeable channel using exact solutions on left and right
% boundaries from Bernardi et al 2023.

% Adapted from main_channel by Kristin Kurianski and Michaela Kubacki Fall
% 2025. Original code: main_channel for Poiseuille flow velocity in a 2D
% channel was developed by Shilpa Khatri and Ricardo Cortez July 2024

% Recent Updates: 
%   - June 2026: Trying different formulations where we use the pressure in
%     step 2 of the channel process.
%   - 3/13: Updated to allow for comparisons with the Segment method
%   - 2/4: Adjusted this file to allow for just having permeability in a
%   middle
%     section of the top and bottom
%   - 12/2: updated all function files to consistently rescale the 
%     Stokeslet and source double matrices by pi
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
N = 160; % Number of source points (along top  and bottom boundaries)
Nx1 = 40*5; % Number of target points in x direction for full channel calc
Nx2 = 80; % Number of target points in y direction for full channel calc

% Setting up channel

%            1 - top (middle portion is perm)
%         +""""""""""""""""""""""""""""""+
%         |                              |
%  3-left |                              | 4-right
% (known) |                              | (known)
%         |                              |                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     |   
%         +""""""""""""""""""""""""""""""+
%         2 - bottom (middle portion is perm)

% Channel Geometry
Lx = pi;
Ly = 1;
xmin = -Lx;
xmax = Lx;
ymin = -Ly;
ymax = Ly;
perm_factor=1;  % will make permeable section  -factor*Lx < x < factor*Lx
scalefactor = 0.05; % for the boundary velocity figures

% Discretization step
ds_x = (xmax - xmin)/N;
ds_y = (ymax - ymin)/( ceil((ymax - ymin)/ds_x));


% Regularization parameter
ep = 5*ds_y;
C = 0.05;
jump_dy = ds_y; % distance from boundary for calculating pressure jump
beta_value = Da*ep*C;
%beta_value = @(x) -0.6 - 0.35*(((pi-x).*(pi+x))/pi^2).^(4/6);
% Discretize Channel Walls (using midpoint to avoid the corners)
% top and bottom
stb = xmin+ds_x/2:ds_x:xmax-ds_x/2;
stb = stb';
% left and right
slr = ymin+ds_y/2:ds_y:ymax-ds_y/2;
slr = slr';

% Define coordinates on each wall (x-coord,y-coord) (source points)
y1_top = stb;   y2_top = ymax*ones(size(stb)); % top wall coordinates (y1_top,y2_top)
y1_bot = stb;   y2_bot = ymin*ones(size(stb)); % bottom wall coordinates (y1_bot,y2_bot)
y1_left = xmin*ones(size(slr));   y2_left = slr; % left wall coordinates (y1_left,y2_left)
y1_right = xmax*ones(size(slr));   y2_right = slr; % right wall coordinates (y1_right,y2_right)
y1 = [y1_top; y1_bot; y1_left; y1_right]; %x-coordinates of all boundary points
y2 = [y2_top; y2_bot; y2_left; y2_right]; %y-coordinates of all boundary points

% Points where velocity will be computed within the channel (target points)
xx1 = linspace(xmin,xmax,Nx1);
xx2 = linspace(ymin,ymax,Nx2);
[x1m, x2m] = ndgrid(xx1, xx2);
x1 = x1m(:); %x-coords of points where computing velocity
x2 = x2m(:); %y-coords of points where computing velocity
dx_g = (xx1(end)-xx1(1))/Nx1;
dy_g = (xx2(end)-xx2(1))/Nx2;

% Indices of permeable region (where the boundary velocity is unknown)
idx = find(-perm_factor*Lx<y1 & y1<perm_factor*Lx);
idx_top = find(-perm_factor*Lx<y1 & y1<perm_factor*Lx & y2==1);
idx_bot = find(-perm_factor*Lx<y1 & y1<perm_factor*Lx & y2==-1);
solidI = find(y1 == xmin | y1 == xmax);
beta = zeros(length(y1),1);
beta(idx) = beta_value;

%%

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



% get exact solution everywhere
[x1gg, x2gg] = meshgrid(xx1, xx2);
[uexact,vexact,~] = permeablechannelexact(x1gg,x2gg,Da);
speed_exact = sqrt(uexact.^2 + vexact.^2);

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

%% Stokeslets Only
%
% Use Stokeslets and exact velocities on top and bottom to find f on Gamma
% Given the exact velocities on boundary of channel (from the true
% solution), find the forces on boundary of channel, then use the boundary
% forces to find the velocities throughout the channel

% All four sides of channel
ftemp = RegStokeslets2D_velocitytoforce([y1,y2], [y1,y2], [u1_bd_exact,u2_bd_exact], ep, mu, blob_num,wt);
% 
% % From boundary forces, calculate channel velocities
ug = RegStokeslets2D_forcetovelocity([y1,y2],ftemp,[x1,x2],ep,mu,blob_num,wt);
% 
ug1 = ug(:,1);
ug2 = ug(:,2);
u1m = reshape(ug1,size(xx1,2),size(xx2,2)); %x-coords of computed velocities
u2m = reshape(ug2,size(xx1,2),size(xx2,2)); %y-coords of computed velocities
% 
% % Errors
u_error = uexact - u1m';
v_error = vexact - u2m';
% 
[MRS_Stokeslet_errors] = calculate_errors(u_error,v_error,xx1,xx2,dx_g,dy_g, 'MRS Stokeslets Only')
% 
plot_streamline(u1m,u2m,x1gg,x2gg,uexact,vexact,y1,y2,'MRS - Stokeslets Only')

% SOMETHING TO POTENTIALLY TRY:

%[p_st] = RegStokeslets2D_forcetovelocity_pressure([y1,y2],ftemp,[x1,x2],ep,mu,blob_num,normals,wt);

% Plot pressure, p_st over channel

% Plot the true pressure against Stokeslet Pressure within
% channel. If matches, try calculating Stokeslet Pressure outside channel.

%% Permeable Channel Problem
% Now we utilize the Stokeslets + Source Doublets approach to simulate the
% permeable channel flow

% 4-Step Permeable Channel Process: Input velocities are the known (exact)
% velocities on left and right, and zero for unknown velocities on top and
% bottom (where it's permeable).

u1 = u1_bd_exact; %x-coordinates of all boundary velocities
u2 = u2_bd_exact; %y-coordinates of all boundary velocities
u1(idx) = 0; 
u2(idx) = 0;

%% STEP 1: Find g-force distribution due to St+SD on solid walls and SD (no
% Stokeslets) on permeable walls

%[g] = RegStokeslets2D_velocityto_gforce_permeable([y1,y2],[y1,y2],...
%    [u1,u2], ep, mu, blob_num, idx, beta, normals, wt);
[g] = RegStokeslets2D_velocityto_gforce_permeable([y1,y2],[y1,y2],...
    [u1,u2], ep, mu, blob_num, idx, beta, normals, wt);
g1 = g(:,1);
g2 = g(:,2);
scaleFactor = 1.5;
figure;
hold on
sk = 1;
quiver(y1(1:sk:end),y2(1:sk:end),g1(1:sk:end),g2(1:sk:end),scaleFactor,'b')
quiver(y1(1:sk:end),y2(1:sk:end),u1(1:sk:end),u2(1:sk:end),scaleFactor,'r')
plot(y1(1:sk:end),y2(1:sk:end),'ro')
axis equal
title('g vs bd velocity')
axis([-4,4,-1.5,1.5])

%% STEP 2A: 
% Use g on full boundary to find pressure on the walls using the Stokeslet
% formula. Then recover the missing velocity on the permeable boundary
% using u dot n = Da [p].

% Calculate Pressure

% Boundary points for permeable section
y1b=y1(idx);
y2b=y2(idx);

% OLD STEP 2 (use source doublets to compute missing boundary velocity)
%[u_perm] = RegStokeslets2D_gtovelocity([y1,y2], g, [y1b,y2b], ep, mu, blob_num, beta, normals, wt);
%u_perm1 = u_perm(:,1);
%u_perm2 = u_perm(:,2);

% Exploration: 
% Calculate Pressure inside and outside channel (wide vertical range). Then
% plot the pressure so we can view the pressure jump.

xxw1 = xx1;
xxw2 = linspace(ymin-1,ymax+1,Nx2); 
[x1wm, x2wm] = ndgrid(xxw1, xxw2);
x1w = x1wm(:); %x-coords of points where computing velocity
x2w = x2wm(:); %y-coords of points where computing velocity

[p] = RegStokeslets2D_gtopressure([y1b,y2b], g, [x1w,x2w], ep, mu, blob_num, beta, normals, wt);
pg = reshape(p,size(xx1,2),size(xxw2,2));

% View a sample vertical slice of the pressure at the 50th x-point
% Can see how it jumps around the top and bottom boundaries (y=+/-1)
figure;
plot(xxw2',pg(50,:)','*-r')
xlabel('y'), ylabel('pressure')
title('Pressure Jump (for fixed x-value)')

% Surface plot pressure (inside and outside channel)
figure; 
surf(x1wm,x2wm,pg)
xlabel('x'),ylabel('y'),zlabel('pressure')
title('Pressure Inside and Outside Channel')

% NEW STEP 2:
% Calculate Pressure outside/inside permeable boundary (top and bottom)
% Set jump_dy under parameters section.

y2b_outside = [y2_top+jump_dy; y2_bot-jump_dy];
[p_outside] = RegStokeslets2D_gtopressure([y1b,y2b], g, [y1b,y2b_outside], ep, mu, blob_num, beta, normals, wt);

y2b_inside = [y2_top-jump_dy; y2_bot+jump_dy];
[p_inside] = RegStokeslets2D_gtopressure([y1b,y2b], g, [y1b,y2b_inside], ep, mu, blob_num, beta, normals, wt);

% Calculate pressure jump
p_jump = p_inside - p_outside;

% Recover missing boundary velocity using
% u dot n = Da * pressure jump and u dot tau = 0
% Note: Can't really justify that we need to divide by 10...

u1(idx) = Da*p_jump .* normals_tb(:,1)/10; 
u2(idx) = Da*p_jump .* normals_tb(:,2)/10;

% Quiver plot of recomputed velocities on boundaries
sk = 1;
figure
subplot(2,1,1);
hold on;
plot(y1,y2,'k.')
hold on
quiver(y1(1:sk:end),y2(1:sk:end),u1(1:sk:end)*scalefactor,u2(1:sk:end)*scalefactor,'b','AutoScale','off')
%quiver(y1(1:sk:end),y2(1:sk:end),u1_bd_exact(1:sk:end)*scalefactor, u2_bd_exact(1:sk:end)*scalefactor,'r','Autoscale','off')
%quiver([y1_left; y1_right],[y2_left; y2_right],u1(203:end),u2(203:end),'b')
title('Computed Velocities on Boundaries after using g')
set(gca,'Fontsize',14)
hold off;
axis equal;
axis([-Lx-1,Lx+1,-1.5,1.5])
subplot(2,1,2);
plot(y1,y2,'k.')
hold on
quiver(y1(1:sk:end),y2(1:sk:end),u1_bd_exact(1:sk:end)*scalefactor, u2_bd_exact(1:sk:end)*scalefactor,'r','Autoscale','off')
hold off
title('Exact Boundary Velocities')
set(gca,'Fontsize',14)
axis equal
axis([-Lx-1,Lx+1,-1.5,1.5])

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
speed = sqrt(u1m.^2 + u2m.^2);

% Errors
u_error = uexact - u1m';
v_error = vexact - u2m';

[MRS_Permeable_errors_pressure] = calculate_errors(u_error,v_error,xx1,xx2,dx_g,dy_g, 'MRS - Permeable')

plot_streamline(u1m,u2m,x1gg,x2gg,uexact,vexact,y1,y2,'MRS - Permeable')

%% Step 2B:  Using pressure jump formulation in integral
% Here we incorporate the pressure jump in the integral formula for the
% Stokeslet pressure solution.

% Calculate Pressure inside and outside channel (wide vertical range)
xxw2 = linspace(ymin-1,ymax+1,Nx2); 
[y1bm, x2wm] = ndgrid(y1b, xxw2);
y1bw = y1bm(:);
x2w = x2wm(:); %y-coords of points where computing pressure

% Calculate Pressure
[p,pjump] =  RegStokeslets2D_gtopressure_integral([y1b,y2b], g, [y1bw,x2w], ep, mu, blob_num, beta, normals, wt);
pg = reshape(p,size(y1b,1),size(xxw2,2));
pjumpg = reshape(pjump,size(y1b,1),size(xxw2,2));

% Plot pressure
figure; 
surf(y1bm,x2wm,pg)
xlabel('x'),ylabel('y'),zlabel('pressure')
title('Pressure Inside and Outside Channel with Jump')

figure; 
surf(y1bm,x2wm,pg-pjumpg)
xlabel('x'),ylabel('y'),zlabel('pressure')
title('Pressure No Jump')

% View a sample vertical slice of the pressure at the 50th x-point
% Can see how it jumps around the top and bottom boundaries (y=+/-1)
figure; hold on
plot(xxw2',pg(1,:)','*-')
plot(xxw2',pg(2,:)','*-')
plot(xxw2',pg(3,:)','*-')
xlabel('y'), ylabel('pressure')
title('Pressure with Jump (for fixed x-values)')

figure; hold on
plot(xxw2',pg(1,:)'-pjumpg(1,:)','*-')
plot(xxw2',pg(2,:)'-pjumpg(2,:)','*-')
plot(xxw2',pg(3,:)'-pjumpg(3,:)','*-')
xlabel('y'), ylabel('pressure')
title('Pressure without Jump (for fixed x-values)')

% Need to find the target points closest to top and bottom boundaries in
% order to calculate the pressure jump. 
% NOTE: Not sure if this is the correct next step here. 

% Find indices in xxw2 closest to the top (y=1) and bottom (y=-1) boundaries

[~, close_top] = mink(abs(xxw2 - 1),2); % closest to y = 1
[~, close_bot] = mink(abs(xxw2 + 1),2); % clostest to y = -1

idx_top_out = max(close_top);
idx_top_in = min(close_top);

idx_bot_out = min(close_bot);
idx_bot_in = max(close_bot);

% Calculate pressure jump

p_jump_top = pg(idx_top,idx_top_in) - pg(idx_top,idx_top_out);
p_jump_bot = pg(idx_bot,idx_bot_in) - pg(idx_bot,idx_bot_out);

% Calculate u using u dot n = Da [p]
% u = (u dot n) n + (u dot tau) tau = (u dot n) n (since u dot tau = 0)
% Again, no justification for why dividing by 10 here makes it a little
% better...
u1(idx_top) = Da*p_jump_top.*normals_top(:,1)/10;
u2(idx_top) = Da*p_jump_top.*normals_top(:,2)/10;

u1(idx_bot) = Da*p_jump_bot.*normals_bot(:,1)/10;
u2(idx_bot) = Da*p_jump_bot.*normals_bot(:,2)/10;

u1(idx) = Da*pjump.*normals_tb(:,1);
u2(idx) = Da*pjump.*normals_tb(:,2);

% Quiver plot of recomputed velocities on boundaries
sk = 1;
figure
subplot(2,1,1);
hold on;
plot(y1,y2,'k.')
hold on
quiver(y1(1:sk:end),y2(1:sk:end),u1(1:sk:end)*scalefactor,u2(1:sk:end)*scalefactor,'b','AutoScale','off')
%quiver(y1(1:sk:end),y2(1:sk:end),u1_bd_exact(1:sk:end)*scalefactor, u2_bd_exact(1:sk:end)*scalefactor,'r','Autoscale','off')
%quiver([y1_left; y1_right],[y2_left; y2_right],u1(203:end),u2(203:end),'b')
title('Computed Velocities on Boundaries after using g')
set(gca,'Fontsize',14)
hold off;
axis equal;
axis([-Lx-1,Lx+1,-1.5,1.5])
subplot(2,1,2);
plot(y1,y2,'k.')
hold on
quiver(y1(1:sk:end),y2(1:sk:end),u1_bd_exact(1:sk:end)*scalefactor, u2_bd_exact(1:sk:end)*scalefactor,'r','Autoscale','off')
hold off
title('Exact Boundary Velocities')
set(gca,'Fontsize',14)
axis equal
axis([-Lx-1,Lx+1,-1.5,1.5])


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
speed = sqrt(u1m.^2 + u2m.^2);

% Errors
u_error = uexact - u1m';
v_error = vexact - u2m';

[MRS_Permeable_errors_pressure_integral] = calculate_errors(u_error,v_error,xx1,xx2,dx_g,dy_g, 'MRS - Permeable')

plot_streamline(u1m,u2m,x1gg,x2gg,uexact,vexact,y1,y2,'MRS - Permeable')


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
%% Streamline Plot function
function plot_streamline(u1m,u2m,x1gg,x2gg,uexact,vexact,y1,y2,method_name)
figure
plot(y1, y2, 'k.')
hold on

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

title(['Computed (blue) and exact (yellow) streamlines - ' method_name])
set(gca, 'FontSize', 20)
xlabel('$x$', 'Interpreter', 'latex')
ylabel('$y$', 'Interpreter', 'latex')
axis equal
axis([-4,4,-1.5,1.5])

end
%% Errors
function [errors] = calculate_errors(u_error,v_error,xx1,xx2,dx,dy, method_name)

u_error_max = max(max(abs(u_error)));
v_error_max = max(max(abs(v_error)));
u_error_L1 = dx*dy*sum(sum(abs(u_error)));
v_error_L1 = dx*dy*sum(sum(abs(v_error)));
u_error_L2 = sqrt(dx*dy*sum(sum(u_error.^2)));
v_error_L2 = sqrt(dx*dy*sum(sum(v_error.^2)));

ColNames = {'Max','L1','L2'};
Emax = [u_error_max; v_error_max];
EL1 = [u_error_L1; v_error_L1];
EL2 = [u_error_L2; v_error_L2];
errors = table(Emax, EL1, EL2,'VariableNames',ColNames);

figure
subplot(2,1,1);
surf(xx1, xx2, u_error, 'EdgeColor', 'none');
colorbar;
title(['Error in u (horizontal) velocity -- ' method_name]);
xlabel('x');
ylabel('y');
axis equal;
view(2); % View from above for 2D representation

subplot(2,1,2);
surf(xx1, xx2, v_error, 'EdgeColor', 'none');
colorbar;
title(['Error in v (vertical) velocity -- ' method_name]);
xlabel('x');
ylabel('y');
axis equal;
view(2); % View from above for 2D representation

end