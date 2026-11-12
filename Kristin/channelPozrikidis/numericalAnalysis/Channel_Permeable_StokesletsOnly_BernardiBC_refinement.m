%Adaptation of Example 3 of Cortez, Fluids 2021
%Channel with inflow and part of membrane permeable on top and bottom

%Developed by Ricardo Cortez, Brittany Leathers, and Michaela Kubacki
%July 2024
%Modified by Kristin Kurianski, Sep 2026

%Find best C1 for eps=Cds^(1/p) for different values of p
%And for different Flow Regimes

clear all 
close all
currDate = datestr(datetime);
currDate = strrep(currDate,'-','_');
currDate = strrep(currDate,' ','_');
currDate = strrep(currDate,':','_');
 
mkdir(currDate)
currPath = [pwd,'/',currDate];
format long


% Set up error matrix and initialize figure
%ErrorMatrix = zeros(length(Nvals),length(epvals));
% Net Flow Error Plot
fig1 = figure;
hold on;
set(gca, 'FontSize', 17);
xlabel('$\epsilon$', 'FontSize', 18,'Interpreter','latex')
ylabel('Error', 'FontSize', 17)
title('Horizontal Velocity Error (L2), FB BC, St only, blob $\psi$',...
    'FontSize', 18,'Interpreter','latex')
hold off;


fig2 = figure;
hold on;
set(gca, 'FontSize', 17);
xlabel('$\epsilon$', 'FontSize', 18,'Interpreter','latex')
ylabel('Error', 'FontSize', 17)
title('Horizontal Velocity Error (Linf), FB BC, St only, blob $\psi$',...
    'FontSize', 18,'Interpreter','latex') 
hold off;


fig3 = figure;
hold on;
set(gca, 'FontSize', 17);
xlabel('$\epsilon$', 'FontSize', 18,'Interpreter','latex')
ylabel('Error', 'FontSize', 17)
title('Vertical Velocity Error (L2), FB BC, St only, blob $\psi$',...
    'FontSize', 18,'Interpreter','latex')
hold off;

fig4 = figure;
hold on;
set(gca, 'FontSize', 17);
xlabel('$\epsilon$', 'FontSize', 18,'Interpreter','latex')
ylabel('Error', 'FontSize', 17)
title('Vertical Velocity Error (Linf), FB BC, St only, blob $\psi$',...
    'FontSize', 18,'Interpreter','latex') 
hold off;


% set up error plots for ep/ds
fig5 = figure;
hold on;
set(gca, 'FontSize', 17);
xlabel('$\epsilon/\Delta s$', 'FontSize', 18,'Interpreter','latex')
ylabel('Error', 'FontSize', 17)
title('Horizontal Velocity Error (L2), FB BC, St only, blob $\psi$',...
    'FontSize', 18,'Interpreter','latex')
hold off;

fig6 = figure;
hold on;
set(gca, 'FontSize', 17);
xlabel('$\epsilon/\Delta s$', 'FontSize', 18,'Interpreter','latex')
ylabel('Error', 'FontSize', 17)
title('Horizontal Velocity Error (Linf), FB BC, St only, blob $\psi$',...
    'FontSize', 18,'Interpreter','latex') 
hold off;

fig7 = figure;
hold on;
set(gca, 'FontSize', 17);
xlabel('$\epsilon/\Delta s$', 'FontSize', 18,'Interpreter','latex')
ylabel('Error', 'FontSize', 17)
title('Vertical Velocity Error (L2), FB BC, St only, blob $\psi$',...
    'FontSize', 18,'Interpreter','latex')
hold off;

fig8 = figure;
hold on;
set(gca, 'FontSize', 17);
xlabel('$\epsilon/\Delta s$', 'FontSize', 18,'Interpreter','latex')
ylabel('Error', 'FontSize', 17)
title('Vertical Velocity Error (Linf), FB BC, St only, blob $\psi$',...
    'FontSize', 18,'Interpreter','latex') 
hold off;

fig9 = figure;
hold on
xlabel('x')
ylabel('velocity')
hold off

%% Parameters to set

%number of points on boundary where velocity is set and force is computed
Nvals=[40 80 160 320 640];

%set Darcy number
Da = 0.4;
%viscosity
mu = 1;
%choose blob
blob_num = 2;
% scale factor for epsilon*ds_y
%ep_factor = 1; % currently using 0.1sqrt(ds)

% Poiseuille flow strength
G = 4;
% chi in PC paper: "dimensionless coeff. determining exit flow rate"
chi = 1;

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

% Necessary for FB solution
% Use Newton's method to compute lambda
lambda = sqrt(2*Da);
fun = @(L) (Da-H/2)*L*cot(L*H)*cos(L*H) + cos(L*H)/2 - H*L*sin(L*H)/2;
lambda = fzero(fun, lambda);

%% Refinement Study

for i=1:length(Nvals)
    N = Nvals(i);

    % Discretization step
    ds_x = (xmax - xmin)/N;
    ds_y = (ymax - ymin)/( ceil((ymax - ymin)/ds_x));
    % Define blob size based on wall discretization
    epvals = logspace(-4,0,20);
    epvals(end-3:end) = [];

    % initialize vectors to store errors
    uL2error = zeros(size(epvals));
    vL2error = zeros(size(epvals));
    uLinferror = zeros(size(epvals));
    vLinferror = zeros(size(epvals));

    for k = 1:length(epvals)
        ep = epvals(k);
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
        u1_left_exact = pois_fun(y2_left);
        u1_right_exact = chi*pois_fun(y2_left);
        u2_left_exact = zeros(size(slr));
        u2_right_exact = zeros(size(slr));

        % Analytical (approx.) velocity in permeable region
        idx_perm = find(c < y1_top & y1_top < L+c);
        u2_top_exact(idx_perm) = -Da*G*sinh(lambda*(y1_top(idx_perm) - xmid))/(lambda*cosh(lambda*L/2));
        u2_bot_exact(idx_perm) = Da*G*sinh(lambda*(y1_top(idx_perm) - xmid))/(lambda*cosh(lambda*L/2));

        % Set zero velocity contribution from stokeslets on permeable parts
        u1_bd_exact = [u1_top_exact; u1_bot_exact; u1_left_exact; u1_right_exact];
        u2_bd_exact = [u2_top_exact; u2_bot_exact; u2_left_exact; u2_right_exact];

        %computing the force
        f = RegStokeslets2D_velocitytoforce([y1,y2], [y1,y2], [u1_bd_exact, u2_bd_exact], ep, mu,...
            blob_num, wt);

        % Calculus velocities across top boundary
        Utop = RegStokeslets2D_forcetovelocity([y1,y2], f, [y1_top, y2_top], ep, mu, blob_num, wt);
        u1_top_comp = Utop(:,1); u2_top_comp = Utop(:,2);

        diffu1 = u1_top_comp - u1_top_exact;
        diffu2 = u2_top_comp - u2_top_exact;
        uL2error(k) = sqrt(ds_x^2*sum(diffu1.^2));
        vL2error(k) = sqrt(ds_y^2*sum(diffu2.^2));
        uLinferror(k) = max(max(abs(diffu1)));
        vLinferror(k) = max(max(abs(diffu2)));

        figure(fig9); hold on;
        quiver(y1_top, y2_top, u1_top_exact,u2_top_exact,'AutoScale', 'off')
        quiver(y1_top, y2_top, u1_top_comp,u2_top_comp,'AutoScale', 'off')
        hold off;

    end

    % Save Current Error results and update figure
    %ErrorMatrix(k,:) = abs(error);
    name = ['N = ' num2str(N)];

    figure(fig1); hold on;
    loglog(epvals, uL2error,'o-', 'LineWidth', 2.5, 'MarkerSize', ...
        10,'DisplayName',name)
    hold off;

    figure(fig2); hold on;
    loglog(epvals,uLinferror,'x--','Linewidth',2.5, 'MarkerSize',10,'DisplayName',name)
    hold off;

    figure(fig3); hold on;
    loglog(epvals, vL2error,'o-', 'LineWidth', 2.5, 'MarkerSize', ...
        10,'DisplayName',name)
    hold off;

    figure(fig4); hold on;
    loglog(epvals,vLinferror,'x--','Linewidth',2.5, 'MarkerSize',10,'DisplayName',name)
    hold off;

    % Save Current Error results and update figure
    % with x-axis ep/ds
    name = ['N = ' num2str(N)];

    figure(fig5); hold on;
    semilogy(epvals/ds_x, uL2error,'o-', 'LineWidth', 2.5, 'MarkerSize', ...
        10,'DisplayName',name)
    hold off;

    figure(fig6); hold on;
    semilogy(epvals/ds_x,uLinferror,'x--','Linewidth',2.5, 'MarkerSize',10,'DisplayName',name)
    hold off;

    figure(fig7); hold on;
    semilogy(epvals/ds_y, vL2error,'o-', 'LineWidth', 2.5, 'MarkerSize', ...
        10,'DisplayName',name)
    hold off;

    figure(fig8); hold on;
    semilogy(epvals/ds_y, vLinferror,'x--','Linewidth',2.5, 'MarkerSize',10,'DisplayName',name)
    hold off;


end

figure(fig1)
set(gca,'YScale','log');
set(gca,'xscale','log');
legend
fig_title = ['Horizontal Error (L2) - no permeability','.fig'];
savefig([currPath,'/',fig_title])

figure(fig2)
set(gca,'YScale','log');
set(gca,'xscale','log');
legend
fig_title = ['Horizontal Error (Linf)- no permeability','.fig'];
savefig([currPath,'/',fig_title])

figure(fig3)
set(gca,'YScale','log');
set(gca,'xscale','log');
legend
fig_title = ['Vertical Error (L2) - no permeability','.fig'];
savefig([currPath,'/',fig_title])

figure(fig4)
set(gca,'YScale','log');
set(gca,'xscale','log');
legend
fig_title = ['Vertical Error (Linf)- no permeability','.fig'];
savefig([currPath,'/',fig_title])

figure(fig5)
set(gca,'YScale','log');
legend
fig_title = ['Horizontal Error (L2) - no permeability - epds','.fig'];
savefig([currPath,'/',fig_title])

figure(fig6)
set(gca,'YScale','log');
legend
fig_title = ['Horizontal Error (Linf)- no permeability - epds','.fig'];
savefig([currPath,'/',fig_title])

figure(fig7)
set(gca,'YScale','log');
legend
fig_title = ['Vertical Error (L2) - no permeability- epds','.fig'];
savefig([currPath,'/',fig_title])

figure(fig8)
set(gca,'YScale','log');
legend
fig_title = ['Vertical Error (Linf)- no permeability - epds','.fig'];
savefig([currPath,'/',fig_title])