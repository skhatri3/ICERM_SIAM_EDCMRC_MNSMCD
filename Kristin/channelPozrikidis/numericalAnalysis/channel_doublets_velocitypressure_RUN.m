% Run MRS with SD for given number of source points around channel N and for
% chosen targert points x1, x2

clear; close all;

% Number of source points
N = 160;

% Input parameters
Da = 0.3; % Darcy number

ep_factor = 0.5;  % ep = ep_factor * ds_y
ep_d_factor = 3; % ep_d = ep_d_factor * ep
beta_factor = 0.34; % beta_value = beta_factor * Da * ep

% Geometry params (pass into function as data structure)
params.L = 3;
params.H = 1;
params.c = 2*params.H;
params.xmin = 0;
params.xmax = params.xmin + params.L + 2*params.c;
params.ymin = 0;
params.ymax = params.ymin + 2*params.H;
params.G = 10;

% Define target points
Nx2 = 80; % Number of target points in y direction for full channel calc
Nx1 = floor(2*pi*Nx2); % Number of target points in x direction for full channel calc
% Points where velocity will be computed within the channel (target points)
xx1 = linspace(params.xmin, params.xmax, Nx1);
xx2 = linspace(params.ymin, params.ymax, Nx2);
[x1m, x2m] = ndgrid(xx1, xx2);
x1 = x1m(:); %x-coords of points where computing velocity
x2 = x2m(:); %y-coords of points where computing velocity

% Compute velocity and pressure
[U, P_sd, P_gauss, y] = channel_doublets_velocitypressure_compute(N, ep_factor, ep_d_factor, Da, beta_factor, params, [x1, x2]);

% Reshape output variables
ug1 = U(:,1);
ug2 = U(:,2);
u1m = reshape(ug1,size(xx1,2),size(xx2,2)); % x-coords of computed velocities everywhere
u2m = reshape(ug2,size(xx1,2),size(xx2,2)); % y-coords of computed velocities everywhere
ummag = sqrt(u1m.^2 + u2m.^2);
pm_sd = reshape(P_sd, size(xx1,2), size(xx2,2));
pm_gauss = reshape(P_gauss, size(xx1,2), size(xx2,2));

% % % Plotting
addpath('..'); 
% method for plots
method_for_plots = "St + SD";

% Create grid for plotting
[x1gg, x2gg] = meshgrid(xx1, xx2);

%Stremlines
figure
plot_streamline_computed_fun(u1m,u2m,x1gg,x2gg,y(:,1), y(:,2),method_for_plots)