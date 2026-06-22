clear all
close all

%% Parameters to set

% Model Parameters
Da = 0.4; % Darcy Number
mu = 1; % Viscosity

% Channel Geometry
Lx = pi;
Ly = 1;
xmin = -Lx;
xmax = Lx;
ymin = -Ly;
ymax = Ly;
perm_factor=1;  % will make permeable section  -factor*Lx < x < factor*Lx

% Numerical Parameters

% Regularization parameter
blob_num = 1; % blob choice
ep_factor = 20; 
C = 0.04; % constant for beta_value 
%scalefactor = 0.05; % for the boundary velocity figures

% number of source and target points
N = 160; % Number of source points (along top  and bottom boundaries)
Nx1 = 40*5; % Number of target points in x direction for full channel calc
Nx2 = 80; % Number of target points in y direction for full channel calc

%% Discretization of channel walls 

% top and bottom
ds_x = (xmax - xmin)/N;
stb = xmin+ds_x/2:ds_x:xmax-ds_x/2;
stb = stb';

% left and right
ds_y = (ymax - ymin)/( ceil((ymax - ymin)/ds_x));
slr = ymin+ds_y/2:ds_y:ymax-ds_y/2;
slr = slr';

% Define coordinates on each wall (y1, y2)
y1_top = stb;   y2_top = ymax*ones(size(stb)); % top wall coordinates (y1_top,y2_top)
y1_bot = stb;   y2_bot = ymin*ones(size(stb)); % bottom wall coordinates (y1_bot,y2_bot)
y1_left = xmin*ones(size(slr));   y2_left = slr; % left wall coordinates (y1_left,y2_left)
y1_right = xmax*ones(size(slr));   y2_right = slr; % right wall coordinates (y1_right,y2_right)
y1 = [y1_top; y1_bot; y1_left; y1_right]; %x-coordinates of all boundary points
y2 = [y2_top; y2_bot; y2_left; y2_right]; %y-coordinates of all boundary points
y1_tb = [y1_top; y1_bot]; % Top and bottom x coordinates
y2_tb = [y2_top; y2_bot]; % Top and bottom y coordinates

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

% Define weights cooresponding to wall coordinates
wt = [ds_x*ones(size(y1_top)); ds_x*ones(size(y1_bot)); ds_y*ones(size(y1_left)); ds_y*ones(size(y1_right))];

% Define blob size based on wall discretization 
ep = ds_y/ep_factor;

%% Setting up permeable region (where the boundary velocity is unknown)

beta_value = Da*ep*C;
idx = find(-perm_factor*Lx<y1 & y1<perm_factor*Lx);
solidI = find(y1 == xmin | y1 == xmax);
beta = zeros(length(y1),1);
beta(idx) = beta_value; %beta is nonzero where permeable 

%% Points where velocity will be computed within the channel

xx1 = linspace(xmin+ds_x/2,xmax-ds_x/2,Nx1);
xx2 = linspace(ymin+ds_y/2,ymax-ds_y/2,Nx2);
[x1gg, x2gg] = meshgrid(xx1, xx2);

[x1m, x2m] = ndgrid(xx1, xx2);
x1 = x1m(:); %x-coords of points where computing velocity
x2 = x2m(:); %y-coords of points where computing velocity

dx_g = (xx1(end)-xx1(1))/Nx1;
dy_g = (xx2(end)-xx2(1))/Nx2;


%% Exact velocities using true solution (Bernardi, 2023)

% Exact velocities on boundaries
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
[uexact,vexact,~] = permeablechannelexact(x1gg,x2gg,Da);
speed_exact = sqrt(uexact.^2 + vexact.^2);

%% Stokeslets Only
%
% Use Stokeslets and exact velocities on top and bottom to find f on Gamma
% Given the exact velocities on boundary of channel (from the true
% solution), find the forces on boundary of channel, then use the boundary
% forces to find the velocities throughout the channel

% All four sides of channel
ftemp = RegStokeslets2D_velocitytoforce([y1,y2], [y1,y2], [u1_bd_exact,u2_bd_exact], ep, mu, blob_num, wt);

% From boundary forces, calculate channel velocities
ug = RegStokeslets2D_forcetovelocity([y1,y2],ftemp,[x1,x2],ep,mu,blob_num,wt);

ug1 = ug(:,1);
ug2 = ug(:,2);
u1m = reshape(ug1,size(xx1,2),size(xx2,2)); %x-coords of computed velocities
u2m = reshape(ug2,size(xx1,2),size(xx2,2)); %y-coords of computed velocities

% Errors
u_error = uexact - u1m';
v_error = vexact - u2m';

[MRS_Stokeslet_errors] = calculate_errors(u_error,v_error,xx1,xx2,dx_g,dy_g, 'MRS Stokeslets Only')

plot_streamline(u1m,u2m,x1gg,x2gg,uexact,vexact,y1,y2,'MRS - Stokeslets Only')



%% Exact solutions (Bernardi, 2023)
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




