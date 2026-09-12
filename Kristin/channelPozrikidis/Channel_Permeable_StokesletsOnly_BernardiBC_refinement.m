%Adaptation of Example 3 of Cortez, Fluids 2021
%Channel with inflow and part of membrane permeable on top and bottom

%Developed by Ricardo Cortez, Brittany Leathers, and Michaela Kubacki
%July 2024
%Modified by Kristin Kurianski, Sep 2026

%Find best C1 for eps=Cds^(1/p) for different values of p
%And for different Flow Regimes

clear all
% close all

%% Parameters to set

%number of points on boundary where velocity is set and force is computed
Nvals=[10 20 40 80 160 320 640 ];

%set Darcy number
Da = 0.2;
%viscosity
mu = 1;
%choose blob
blob_num = 2;
% scale factor for epsilon*ds_y
%ep_factor = 1; % currently using 0.1sqrt(ds)

% Poiseuille flow strength
G = 2;
% chi in PC paper: "dimensionless coeff. determining exit flow rate"
chi = 1;

% Channel geometry
L = 2; % length of permeable portion of the channel
H = 1; % radius of the channel
c = 2*H; % extension length of the channel

xmin = -L - c;
xmax = L + c;
ymin = -H;
ymax = H;

% needed for FB solution
% Use Newton's method to compute lambda
lambda = sqrt(2*Da);
fun = @(L) (Da-H/2)*L*cot(L*H)*cos(L*H) + cos(L*H)/2 - H*L*sin(L*H)/2;
lambda = fzero(fun, lambda);

eumaxnorm_away=zeros(length(Nvals),1);
eu2norm_away=zeros(length(Nvals), 1);
eu_pointaway=zeros(length(Nvals), 1);
eumaxnorm_near=zeros(length(Nvals),1);
eu2norm_near=zeros(length(Nvals), 1);

%% Refinement Study

for i=1:length(Nvals)
    N = Nvals(i)

    % Discretization step
    ds_x = (xmax - xmin)/N;
    ds_y = (ymax - ymin)/( ceil((ymax - ymin)/ds_x));

    % Define blob size based on wall discretization
    %ep = ds_y*ep_factor;
    %ep=0.1*ds_y^(1/2);
    ep = ds_y;


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
    
    % % %unit normals
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

    %Quadrature weights
    % Currently using midpoint rule.
    wt = [ds_x*ones(size(y1_top)); ds_x*ones(size(y1_bot)); ds_y*ones(size(y1_left)); ds_y*ones(size(y1_right))];

    % No slip at top and bottom
    % (set to 0, although technically unknown in permeable region)
    u1_top_exact = zeros(size(y1_top));
    u2_top_exact = zeros(size(y2_top));
    u1_bot_exact = zeros(size(y1_bot));
    u2_bot_exact = zeros(size(y2_bot));
    % Poiseuille flow at inlet and outlet
    u_pois = -G*(y2_left - ymin).*(y2_left - ymax);%/2/mu;
    u1_left_exact = u_pois;
    u1_right_exact = chi*u_pois;
    u2_left_exact = zeros(size(slr));
    u2_right_exact = zeros(size(slr));
    
    % Analytical (approx.) velocity in permeable region
    idx_perm = find(-L < y1_top & y1_top < L);
    u2_top_exact(idx_perm) = -Da*G*sinh(lambda*(y1_top(idx_perm)))/(lambda*cosh(lambda*L));
    u2_bot_exact(idx_perm) = Da*G*sinh(lambda*(y1_top(idx_perm)))/(lambda*cosh(lambda*L));

    % Set zero velocity contribution from stokeslets on permeable parts
    u1_bd_exact = [u1_top_exact; u1_bot_exact; u1_left_exact; u1_right_exact];
    u2_bd_exact = [u2_top_exact; u2_bot_exact; u2_left_exact; u2_right_exact];

    %computing the force
    f = RegStokeslets2D_velocitytoforce([y1,y2], [y1,y2], [u1_bd_exact, u2_bd_exact], ep, mu,...
        blob_num, wt);
    
    %calculate on grid away from the boundary
    %calculate on grid away from the boundary
    xmin_away = -L/2;
    ymin_away = -0.1;
    Lx_away = L; xmax_away = xmin_away + Lx_away;
    Ly_away = 0.2; ymax_away = ymin_away + Ly_away;
    dx_away = (ymax_away-ymin_away)/(2^6);
    % dx=ds;
    % dxs(i)=dx;
    Nx_away = round(Lx_away/dx_away);
    Ny_away = round(Ly_away/dx_away);
    xg_away = dx_away*(0:Nx_away-1) + xmin_away;
    yg_away = dx_away*(0:Ny_away-1) + ymin_away;
    [xg_away,yg_away]=ndgrid(xg_away,yg_away);
    xgv_away=reshape(xg_away, Nx_away*Ny_away,1);
    ygv_away=reshape(yg_away, Nx_away*Ny_away,1);
    
    % Compute velocities on grid
    Ugrid = RegStokeslets2D_forcetovelocity([y1,y2],f,[xgv_away,ygv_away],ep,mu,blob_num,wt);
    ug_away=reshape(Ugrid(:,1), Nx_away, Ny_away);
    vg_away=reshape(Ugrid(:,2), Nx_away, Ny_away);

    % !!NOTE: Use ST sol as "true" and SD solution as approx. Need to edit
    % this!!
    u_soln_away = permeableChannelExact(xg_away, yg_away, Da);
    uerror_away = abs(ug_away-u_soln_away);
    eumaxnorm_away(i) = max(uerror_away(:));
    eu2norm_away(i) = sqrt(dx_away^2*sum(sum(uerror_away.^2)));

end
%%

colorp=[0.4940, 0.1840, 0.5560];
colorlb=[0.3010, 0.6450, 0.9930];
colorg=[0.4660, 0.6740, 0.1880];
colordb=	[0, 0.4470, 0.7410];

forploty0 = 0.1;
forplotdx0 = 0.02;
dxforplot = 1./Nvals;
yforplot = forploty0/forplotdx0^1*dxforplot.^(2);
forploty02 = 25;
forplotdx02 = 10;
dxforplot2 = 1./Nvals;
yforplot2 = forploty02/forplotdx02^2*dxforplot2.^1;

%Refinement study plots
figure;
loglog(Nvals, eumaxnorm_away,'o-', 'LineWidth', 2.5, 'MarkerSize', ...
    10, 'Color', colorlb)
hold on
loglog(Nvals, yforplot, 'LineWidth', 2,'Color', colorg)
text(100,1/7000,'$\Delta s^{2}$', 'Color', colorg, 'FontSize', 14, ...
    'Interpreter', 'latex');
hold off
axis([10 1000 10^(-7) 10^(0) ]);
set(gca, 'FontSize', 17);
xlabel('$N$', 'FontSize', 18,'Interpreter','latex')
ylabel('$||e||_{\infty}$', 'FontSize', 17, 'Interpreter','latex')
title('Refinement Study, $\epsilon=0.1ds^{1/2}$, "away"',...
    'FontSize', 18,'Interpreter','latex')
