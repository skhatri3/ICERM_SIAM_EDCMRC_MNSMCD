% Run MRS with SD for given number of source points around channel N and for
% chosen targert points x1, x2

clear; close all;

% Number of source points
N = 160;

% Input parameters
params.Da = 0.3; % Darcy number

ep_factor = 0.034;  % ep = ep_factor * ds_y

% Geometry params (pass into function as data structure)
params.L = 3;
params.H = 1;
params.c = 2*params.H;
params.xmin = 0;
params.xmax = params.xmin + params.L + 2*params.c;
params.ymin = 0;
params.ymax = params.ymin + 2*params.H;
params.G = 10;

% Use Newton's method to compute eta eigenvalue
params.eta = sqrt(2*params.Da);
fun = @(L) (params.Da-params.H/2)*L*cot(L*params.H)*cos(L*params.H) + cos(L*params.H)/2 - params.H*L*sin(L*params.H)/2;
params.eta = fzero(fun, params.eta);

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
[U, P, y] = channel_stokeslets_velocitypressure_compute(N, ep_factor, params, [x1, x2]);

% Reshape output variables
ug1 = U(:,1);
ug2 = U(:,2);
u1m = reshape(ug1,size(xx1,2),size(xx2,2)); % x-coords of computed velocities everywhere
u2m = reshape(ug2,size(xx1,2),size(xx2,2)); % y-coords of computed velocities everywhere
ummag = sqrt(u1m.^2 + u2m.^2);
pm_mrs = reshape(P, size(xx1,2), size(xx2,2));

%% Error calculations
% Create grid for calculations
%if isempty(x1gg) || isempty(x2gg)
    [x1gg, x2gg] = meshgrid(xx1, xx2);
%end

% Compute exact pressure
% Need iscretization step
ds_x = (params.xmax - params.xmin)/N;
ds_y = (params.ymax - params.ymin)/( ceil((params.ymax - params.ymin)/ds_x));

p_exact = exactPressureFun_Fourier(x1gg, x2gg, params);

[errorTable, errorMat] = calculate_errors_fun([], [], pm_mrs, p_exact', ds_x, ds_y);

% plot
figure
pc = pcolor(x1gg,x2gg,p_exact);
set(pc,'facealpha',0.8)
shading interp
hold on
plot(y(:,1), y(:, 2), '.k')
axis equal


%% Plotting
% Create grid for calculations
[x1gg, x2gg] = meshgrid(xx1, xx2);

addpath('..');
addpath('../plotFuns'); 

% Streamlines (need grid points, x1gg and x2gg)
figure
plot_streamline_computed_fun(u1m, u2m, x1gg, x2gg, y(:,1), y(:,2), 'Stokeslets Only')
hold on
% indicate boundaries of permeable region
plot([params.c, params.c], [0, 2*params.H], '--k', [params.L + params.c, params.L + params.c], [0, 2*params.H], '--k', 'LineWidth', 1.3)
text(params.c - 0.1, -0.15, '$c$', 'Interpreter','latex', 'FontSize', 16)
text(params.L + params.c - 0.3, -0.2, '$L+c$', 'Interpreter','latex', 'FontSize', 16)
% formatting
ax = gca;
ax.FontSize = 16;


%% Error throughout channel
figure
plot(y(:,1),y(:,2),'k.')
hold on
abs_pressure_error = abs(p_exact - pm_mrs');
pc = pcolor(x1gg, x2gg, abs_pressure_error);
set(pc,'facealpha',0.8)
shading interp
colorbar('southoutside')
set(gca,'FontSize',18)
axis equal
title('Pressure error throughout channel - Stokeslets Only - Absolute Error')