% Use bisection method to find epsilon that gives correct flux across half
% upper boundary
clear
close all

addpath('./stokeslets_codes');

%set Darcy number
Da = 0.4;
%viscosity
mu = 1;
%choose blob
blob_num = 2;

% Poiseuille flow strength
G = 4;
% chi in PC paper: "dimensionless coeff. determining exit flow rate"
chi = 1;

% Number of source and target points
N = 160; % Number of source points (along top  and bottom boundaries)
Nx2 = 40; % Number of target points in y direction for full channel calc
Nx1 = floor(pi*Nx2); % Number of target points in x direction for full channel calc

% Setting forces and computing velocity
% Channel geometry
L = 3; % length of permeable portion of the channel
H = 1; % radius of the channel
c = 2*H; % extension length of the channel
xmin = 0;
xmax = L + 2*c;
xmid = (L+2*c)/2;
ymin = 0;
ymax = 2*H;

% Poiseuille flow function for inlet and outlet
pois_fun = @(y) -G*(y - ymin).*(y - ymax)/2;

% Discretization step
ds_x = (xmax - xmin)/N;
ds_y = (ymax - ymin)/( ceil((ymax - ymin)/ds_x));

% Define initial blob size based on wall discretization
ep_min = ds_y/100;
ep_max = 10*ds_y;
ep0 = (ep_min + ep_max)/2;
ep_vec = ep0;

% Necessary for FB solution
% Use Newton's method to compute lambda
lambda = sqrt(2*Da);
fun = @(L) (Da-H/2)*L*cot(L*H)*cos(L*H) + cos(L*H)/2 - H*L*sin(L*H)/2;
lambda = fzero(fun, lambda);
% define the functions g(y) and g'(y) from the notes
g=@(y) (Da - H/2)*cot(lam*H)*sin(lam*y) + y/2.*cos(lam*y);
gp=@(y) (Da - H/2)*lam*cot(lam*H)*cos(lam*y) + 1/2.*cos(lam*y) - lam*y/2.*sin(lam*y);

% Compute theoretical flux across half upper boundary
flux_theory = -Da*G*(1-cosh(lambda*L/2))/(lambda^2*cosh(lambda*L/2));
flux_numerical = 2*flux_theory;

% Tolerance for bisection method
tol = 1e-2;

while abs(flux_numerical - flux_theory) > tol

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

    % Points where velocity will be computed within the channel (target points)
    xx1 = linspace(xmin,xmax,Nx1);
    xx2 = linspace(ymin,ymax,Nx2);
    [x1m, x2m] = ndgrid(xx1, xx2);
    x1 = x1m(:); %x-coords of points where computing velocity
    x2 = x2m(:); %y-coords of points where computing velocity
    % Define quadrature weights cooresponding to wall coordinates.
    % Currently using midpoint rule.
    % Weights for targets/sources at corners and just outside inlet/outlet are
    % set to ds_y (possibly not correct) but will be divided out so shouldn't cause issues.
    wt = [ds_x*ones(size(y1_top)); ds_x*ones(size(y1_bot)); ds_y*ones(size(y1_left)); ds_y*ones(size(y1_right))];

    % No slip at top and bottom
    % (set to 0, although technically unknown in permeable region)
    u1_top_exact = zeros(size(y1_top));
    u2_top_exact = zeros(size(y2_top));
    u1_bot_exact = zeros(size(y1_bot));
    u2_bot_exact = zeros(size(y2_bot));
    % Poiseuille flow at inlet and outlet
    u1_left_exact = pois_fun(y2_left);
    u1_right_exact = chi*pois_fun(y2_left);
    u2_left_exact = zeros(size(slr));
    u2_right_exact = zeros(size(slr));

    % Analytical (approx.) velocity in permeable region
    idx_perm = find(c < y1_top & y1_top < L+c);
    u2_top_exact(idx_perm) = -Da*G*sinh(lambda*(y1_top(idx_perm)-xmid))/(lambda*cosh(lambda*L/2));
    u2_bot_exact(idx_perm) = Da*G*sinh(lambda*(y1_top(idx_perm)-xmid))/(lambda*cosh(lambda*L/2));

    % Set zero velocity contribution from stokeslets on permeable parts
    u1_bd_exact = [u1_top_exact; u1_bot_exact; u1_left_exact; u1_right_exact];
    u2_bd_exact = [u2_top_exact; u2_bot_exact; u2_left_exact; u2_right_exact];

    % Compute forces
    f_temp = RegStokeslets2D_velocitytoforce_KK([y1,y2], [y1,y2], [u1_bd_exact,u2_bd_exact], ep0, mu, blob_num,wt);

    % % Compute velcoity using Bernardi BCs % %
    ug = RegStokeslets2D_forcetovelocity([y1,y2],f_temp,[x1,x2],ep0,mu,blob_num,wt);
    ug1 = ug(:,1);
    ug2 = ug(:,2);
    u1m = reshape(ug1,size(xx1,2),size(xx2,2)); % x-coords of computed velocities everywhere
    u2m = reshape(ug2,size(xx1,2),size(xx2,2)); % y-coords of computed velocities everywhere
    ummag = sqrt(u1m.^2 + u2m.^2);

    % Numerical flux across half the top boundary
    idx_flux = find(y1_top < c + L/2);
    y1_flux = y1_top(idx_flux);
    normals_top = zeros(length(y1_flux),2);
    normals_top(:,2) = 1;
    flux_numerical = ds_x*sum(dot(normals_top, [zeros(size(y1_flux)), u2_top_exact(1:length(y1_flux))]));

    if flux_theory - flux_numerical < -tol % too much fluid leaving
        ep_max = ep0;
        ep0 = (ep_min + ep0)/2;
        ep_vec = [ep_vec, ep0];
    elseif flux_theory - flux_numerical > tol % not enough fluid leaving
        ep_min = ep0;
        ep0 = (ep_max + ep0)/2;
        ep_vec = [ep_vec, ep0];
    end

end