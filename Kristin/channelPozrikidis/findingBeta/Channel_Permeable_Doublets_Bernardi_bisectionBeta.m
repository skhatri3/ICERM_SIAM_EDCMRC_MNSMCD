% Pozrikidis geometry, Berndardi boundary condition
% Determine beta given Da, epsilon, and N

clear all 
close all

addpath('../stokeslets_codes'); 

%% Parameters
mu = 1; % viscosity
blob_num = 2; %blob choice
G = 4; % poiseuille flow strength
Da = 0.3; % Darcy number
chi = 1; % scale factor for exit velocity

%channel geometry
L = 3; % length of permeable region
H = 1; % height of channel
c = 2*H; % length of nonpermeable region on one side (scale factor of H based on Pozrikidis)
xmin = 0;
xmax = xmin + L + 2*c;
ymin = 0;
ymax = ymin + 2*H;
perm_min = c;         %start of permeable part
perm_max = L + c;  %end of permeable part

% Poiseuille flow function for inlet and outlet
pois_fun = @(y) -G*(y - ymin).*(y - ymax)/2;

%discretization of channel
N = 160; %number of points on boundary where velocity is set and force is computed 
Nx2 = 40; % Number of target points in y direction for full channel calc
Nx1 = floor(pi*Nx2); % Number of target points in x direction for full channel calc
% Discretization step
ds_x = (xmax - xmin)/N;
ds_y = (ymax - ymin)/( ceil((ymax - ymin)/ds_x));

%regularization parameter
ep1 = 0.010;
ep2 =ep1;

% Necessary for FB solution
% Use Newton's method to compute eta
eta = sqrt(2*Da);
fun = @(L) (Da-H/2)*L*cot(L*H)*cos(L*H) + cos(L*H)/2 - H*L*sin(L*H)/2;
eta = fzero(fun, eta);

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

%indices of permeable region
permidx = find(y1 < L + c & y1 > c);
% y's for permeable region
y1b=y1(permidx);
y2b=y2(permidx);

% Target Coordinates for outlet
x2_out= (ymin + ds_y/2:ds_y:ymax - ds_y/2)';
x1_out= xmax*ones(size(x2_out));

% Weights for targets/sources at corners and just outside inlet/outlet are
% set to ds_y (possibly not correct) but will be divided out so shouldn't cause issues.
wt = [ds_x*ones(size(y1_top)); ds_x*ones(size(y1_bot)); ds_y*ones(size(y1_left)); ds_y*ones(size(y1_right))];

%% Bisection
% Compute analytical flux across top half of permeabl portion
betaUpper = 0.05;
betaLower = 0.00009;

Rtarget = Da*G*(cosh(eta*L/2)-1)/eta^2/cosh(eta*L/2); % analytical flux across half of permeable top boundary

%*** !!! need to actually compute these !!! ***
RoutbLower = 0.8;
RoutbUpper = 1.5;

RdiffLower = RoutbLower - Rtarget; % for current ep and betaUpper
RdiffUpper = RoutbUpper - Rtarget; % for current ep and betaLower

tol = 1e-6;
itMax = 30;
it_count = 0;

for k = 1:itMax
    it_count = it_count + 1;
    betaCurrent = (betaLower + betaUpper)/2;

    % Set boundary velocities
    % No slip at top and bottom
    % (set to 0, although technically unknown in permeable region)
    u1_top_exact = zeros(size(y1_top)); u2_top_exact = zeros(size(y2_top));
    u1_bot_exact = zeros(size(y1_bot)); u2_bot_exact = zeros(size(y2_bot));
    % Poiseuille flow at inlet and outlet
    u1_left_exact = pois_fun(y2_left); % call poiseuille flow profile function
    u1_right_exact = chi*u1_left_exact; % scale poiseuille flow at outlet
    u2_left_exact = zeros(size(slr)); % purely horizontal flow at inlet
    u2_right_exact = zeros(size(slr)); % purely horizontal flow at outlet

    % Set zero velocity contribution from stokeslets on permeable parts
    u1 = [u1_top_exact; u1_bot_exact; u1_left_exact; u1_right_exact];
    u2 = [u2_top_exact; u2_bot_exact; u2_left_exact; u2_right_exact];

    beta = zeros(length(y1),1);
    beta(permidx) = betaCurrent;
    y1b = y1(permidx);
    y2b = y2(permidx);

    % Find g forces due to known velocities
    [g] = RegStokeslets2D_velocityto_gforce_SK([y1,y2], [y1,y2], [u1,u2],ep1,ep2,mu,blob_num, permidx, beta, normals, wt);
    
    % Find permeable velocity due to g forces
    [u_perm] = RegStokeslets2D_gtovelocity([y1,y2], g, [y1b,y2b], ep2, mu, blob_num, beta, normals, wt);
    u_perm1 = u_perm(:,1);
    u_perm2 = u_perm(:,2);

    % Replace the placeholder zeros in u1,u2 with the velocities
    % [zeros(size(u1_top_exact)); zeros(size(u1_bot_exact)); u1_left_exact; u1_right_exact; zeros(size(u1_corner_exact))]; %x-coordinates of all boundary velocities
    u1(permidx) = u_perm1;
    u2(permidx) = u_perm2;

    % Compute forces on full boundary with new known boundary velocities
    ftemp = RegStokeslets2D_velocitytoforce_KK([y1,y2], [y1,y2], [u1,u2], ep1, mu, blob_num,wt);
    f1 = ftemp(:,1);
    f2 = ftemp(:,2);

    %outlet velocity
    uu = RegStokeslets2D_forcetovelocity([y1,y2],ftemp,[x1_out,x2_out],ep1,mu,blob_num,wt);


    % Numerical flux across half the top boundary
    %idx_flux = find(y1_top >= c & y1_top < c+L/2);
    %y1_flux = y1_top(idx_flux);

    %normals_top=zeros(length(y1_flux),2);
    %normals_top(:,2)=1;

    y1b_tophalf = y1b(1:floor(end/2));
    idx_permflux = find(y1b_tophalf < c+L/2);

    normals_tophalf=zeros(length(y1_flux),2);
    normals_tophalf(:,2)=1;

    Rtop = ds_x*sum(dot(normals_tophalf, [zeros(size(y1_flux)), u_perm2(idx_permflux)]));

    Rdiff = Rtop - Rtarget; % current difference in outflows

    current = [it_count, betaUpper, betaLower, RdiffLower];
    disp(current)

    if (abs(Rdiff) < tol)
        disp('Bisection Converged')
        break
    end

    if sign(RdiffLower) ~= sign(RdiffLower) % root is in [betaLower, betaCurr]
        betaUpper = betaCurrent;
        RdiffUpper = Rdiff;
    else % root is in [betaCurrent, betaUpper]
        betaLower = betaCurrent;
        RdiffLower = Rdiff;
    end

    if (it_count >= itMax)
        disp('Did not converge within tol. Exceeded itMax')
    end

end