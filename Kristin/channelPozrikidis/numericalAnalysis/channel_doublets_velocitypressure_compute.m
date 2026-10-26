% Use MRS+SD to compute u at desired target locations given fluid parameters
function [U, P_st, P_gauss, y] = channel_doublets_velocitypressure_compute(N, ep_factor, ep_d_factor, beta_factor, params, x)

% INPUTS:
% N -- number of source points
% ep_factor -- scale factor for blob size for stokeslets; multiplies ds_y
% ep_d_factor -- scale factor for blob size for doublets; multiplies ep
% Da -- Darcy number
% beta_factor -- scale factor for beta value; multiplies Da*ep
% params -- data structure containing geometry parameters
% x1, x2 -- Target points x-comp and y-comp, resp.

% OUTPUTS:
% U -- Velocity array U with x-component in column 1 and y-component in
%               column 2 at target points
% P_st -- Pressure at target points using stokeslets pressure expression
% P_gauss -- Pressure at target points using gauss formulation
% y -- source points (x-comp in column 1, y-comp in column 2)

addpath('../stokeslets_codes'); 

% Parameters that will be unchanged across runs
%viscosity
mu = 1; 
% Regularization parameter
blob_num = 2; % blob choice

% Poiseuille flow strength
G = params.G; 
% chi in PC paper: "dimensionless coeff. determining exit flow rate"
chi = 1;

% Setting forces and computing velocity 
% Channel geometry
L = params.L; % length of permeable portion of the channel
H = params.H; % radius of the channel
c = params.c; % extension length of the channel
xm = (L + 2*c)/2; % midpoint of channel
Da = params.Da; % Darcy number

% Use Newton's method to compute eta eigenvalue
eta = sqrt(2*Da);
fun = @(L) (Da-H/2)*L*cot(L*H)*cos(L*H) + cos(L*H)/2 - H*L*sin(L*H)/2;
eta = fzero(fun, eta);

xmin = params.xmin;
xmax = params.xmax;
ymin = params.ymin;
ymax = params.ymax;

% Poiseuille flow function for inlet and outlet
pois_fun = @(y) -G*(y - ymin).*(y - ymax)/2;

% Discretization step
ds_x = (xmax - xmin)/Nx1;
ds_y = (ymax - ymin)/( ceil((ymax - ymin)/ds_x));

% Define blob size based on wall discretization 
ep = ds_y*ep_factor;
ep_d =ep*ep_d_factor; 

beta_value = Da*ep*beta_factor;

%discretization step of top/bottom walls
stb = xmin+ds_x/2:ds_x:xmax-ds_x/2;
stb = stb';
%discretization step of left/right wall (inlet/outlet)
slr = ymin+ds_y/2:ds_y:ymax-ds_y/2;
slr = slr';

% Define coordinates on each wall (x-coord,y-coord) (source points)
y1_top = stb;   y2_top = ymax*ones(size(stb)); % top wall coordinates (y1_top,y2_top)
y1_bot = stb;   y2_bot = ymin*ones(size(stb)); % bottom wall coordinates (y1_bot,y2_bot)
y1_left = xmin*ones(size(slr));   y2_left = slr; % left wall coordinates (y1_left,y2_left)
y1_right = xmax*ones(size(slr));   y2_right = slr; % right wall coordinates (y1_right,y2_right)
y1 = [y1_top; y1_bot; y1_left; y1_right]; %x-coordinates of all boundary points
y2 = [y2_top; y2_bot; y2_left; y2_right]; %y-coordinates of all boundary points
y = [y1, y2];

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

% Define quadrature weights cooresponding to wall coordinates.
% Currently using midpoint rule.
% Weights for targets/sources at corners and just outside inlet/outlet are
% set to ds_y (possibly not correct) but will be divided out so shouldn't cause issues.
wt = [ds_x*ones(size(y1_top)); ds_x*ones(size(y1_bot)); ds_y*ones(size(y1_left)); ds_y*ones(size(y1_right))];

% No slip at top and bottom
u1_top_exact = zeros(size(y1_top)); u2_top_exact = zeros(size(y2_top));
u1_bot_exact = zeros(size(y1_bot)); u2_bot_exact = zeros(size(y2_bot));
% Poiseuille flow at inlet and outlet
u1_left_exact = pois_fun(y2_left);
u1_right_exact = chi*u1_left_exact;
u2_left_exact = zeros(size(slr));
u2_right_exact = zeros(size(slr));

% Set zero velocity contribution from stokeslets on permeable parts
u1 = [u1_top_exact; u1_bot_exact; u1_left_exact; u1_right_exact];
u2 = [u2_top_exact; u2_bot_exact; u2_left_exact; u2_right_exact];

idx = find(c < y1 & y1 < L+c);
beta = zeros(length(y1),1);
beta(idx) = beta_value;

% Find g-force distribution due to St+SD on solid walls and SD (no Stokeslets) on permeable walls
[g] = RegStokeslets2D_velocityto_gforce_SK([y1,y2], [y1,y2], [u1,u2],ep,ep_d,mu,blob_num, idx, beta, normals, wt); 

% Boundary points for permeable section
y1b=y1(idx);
y2b=y2(idx);
% Inputs to function are source, g force, target, ...
[u_perm] = RegStokeslets2D_gtovelocity([y1,y2], g, [y1b,y2b], ep_d, mu, blob_num, beta, normals, wt);
u_perm1 = u_perm(:,1);
u_perm2 = u_perm(:,2);

% Replace the placeholder zeros in u1,u2 with the velocities
u1(idx) = u_perm1;
u2(idx) = u_perm2;


% Use boundary velocities to solve for forces on boundary
% Now that the missing velocities are recovered, the next steps are the regular
% "Stokeslet" problem, given boundary velocities, find boundary forces,
% then find channel velocities.
ftemp = RegStokeslets2D_velocitytoforce_KK([y1,y2], [y1,y2], [u1,u2], ep, mu, blob_num,wt);

% Compute velocity everywhere using the computed forces
% Unpack target points
x1 = x(:,1); 
x2 = x(:,2);

% Compute pressure using gauss formulation
% Define tangent vectors
tangent_top = zeros(length(stb),2); % unit tangent for top:
tangent_top(:,1) = -1;
tangent_bot = zeros(length(stb),2); % unit tangent for bottom
tangent_bot(:,1) = 1;
tangent_left = zeros(length(slr),2); % unit tangent for left  wall 
tangent_left(:,2) = -1;
tangent_right = zeros(length(slr),2); % unit tangent for right wall 
tangent_right(:,2) = 1;
% tangent on full boundary
tangents = [tangent_top; tangent_bot; tangent_left; tangent_right];
% Pressure computation
P_gauss = RegStokeslets2D_forcetopressure_gauss([y1,y2], ftemp, [x1, x2], ep, mu, beta, blob_num, normals, tangents, wt);

% Compute velocity everywhere using the computed forces
U = RegStokeslets2D_forcetovelocity([y1,y2],ftemp,[x1,x2],ep,mu,blob_num,wt);

% Computed pressure using SD
P_st = RegStokeslets2D_forcetopressure([y1,y2], ftemp, [x1, x2], ep, mu, blob_num, wt);