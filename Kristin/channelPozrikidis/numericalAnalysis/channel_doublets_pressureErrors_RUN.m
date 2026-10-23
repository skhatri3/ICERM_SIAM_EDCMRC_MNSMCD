% Run MRS for given number of source points around channel N and for
%       chosen targert points x1, x2
% Computes errors in pressure based on exact pressure computation using
%       Fourier series

clear; close all

% Create directory for saving figures
dformat = datetime; % Put date and time in directory name
dformat.Format = "ddMMMyyyy_'time_h'HH'm'mm's'ss";
directoryName = ['PressureError_', char(dformat)];
mkdir(directoryName)
currPath = [pwd,'/',directoryName];
format long

% Number of source points
N_vec = [20, 40, 80, 160, 320];

% Determine which target points to use in error calculation
% Options: 'midline', 'permStart', 'permEnd'
targetType = 'midline';

% initialize vectors to store errors
pL2error = zeros(size(N_vec));
pLinferror = zeros(size(N_vec));

% Input parameters
Da = 0.3; % Darcy number

% Geometry params (pass into function as data structure)
params.L = 3;
params.H = 1;
params.c = 2*params.H;
params.xmin = 0;
params.xmax = params.xmin + params.L + 2*params.c;
params.ymin = 0;
params.ymax = params.ymin + 2*params.H;
params.G = 10;

% Use Newton's method to compute eta eigenvalue for true solution
params.eta = sqrt(2*Da);
fun = @(L) (Da-params.H/2)*L*cot(L*params.H)*cos(L*params.H) + cos(L*params.H)/2 - params.H*L*sin(L*params.H)/2;
params.eta = fzero(fun, params.eta);


% Set up error matrix and initialize figure
%ErrorMatrix = zeros(length(Nvals),length(epvals));
% Net Flow Error Plot
fig1 = figure;
hold on;
set(gca, 'FontSize', 17);
xlabel('$\epsilon$', 'FontSize', 18,'Interpreter','latex')
ylabel('Error', 'FontSize', 17)
title('Pressure Error (L2), Permeable, Stokeslets Only, Bob $\psi$',...
    'FontSize', 18,'Interpreter','latex')
hold off;

fig2 = figure;
hold on;
set(gca, 'FontSize', 17);
xlabel('$\epsilon$', 'FontSize', 18,'Interpreter','latex')
ylabel('Error', 'FontSize', 17)
title('Pressure Error (Linf), Permeable, Stokeslets Only, Blob $\psi$',...
    'FontSize', 18,'Interpreter','latex') 
hold off;

% set up error plots for ep/ds
fig3 = figure;
hold on;
set(gca, 'FontSize', 17);
xlabel('$\epsilon/\Delta s$', 'FontSize', 18,'Interpreter','latex')
ylabel('Error', 'FontSize', 17)
title('Pressure Error (L2), Permeable, Stokeslets Only, Blob $\psi$',...
    'FontSize', 18,'Interpreter','latex')
hold off;

fig4 = figure;
hold on;
set(gca, 'FontSize', 17);
xlabel('$\epsilon/\Delta s$', 'FontSize', 18,'Interpreter','latex')
ylabel('Error', 'FontSize', 17)
title('Pressure Error (Linf), Permeable, Stokeslets Only, Blob $\psi$',...
    'FontSize', 18,'Interpreter','latex') 
hold off;

for j =1:length(N_vec) % loop over number of source points
    % Current number of source points
    N = N_vec(j);

    % % Define target points
    % get discretization for target points based on current N
    ds_y = (params.ymax - params.ymin)/N;
    ds_x = (params.xmax - params.xmin)/( ceil((params.xmax - params.xmin)/ds_y));
    
    % Target point choice based on targetType defined before loops
    if strcmp(targetType, 'midline')
        x2 = params.ymin+ds_y/2:ds_y:params.ymax-ds_y/2;
        x2 = x2(:);
        x1 = ((params.xmin + params.xmax)/2)*ones(size(x2));
    end

    % blob sizes for loop
    epvals = logspace(-4,0,20);

    for k = 1:length(epvals) % loop over blob size values
        % Current epsilon value for blob size
        ep = epvals(k);

        % Compute velocity and pressure at chosen target points
        [U, P, y] = channel_stokeslets_velocitypressure_compute(N, ep, Da, params, [x1, x2]);
        u1_comp = U(:,1);
        u2_comp = U(:,2);
        umag_comp = sqrt(u1_comp.^2 + u2_comp.^2);
        p_comp = P;
        
        % Compute exact pressure
        p_exact = exactPressureFun(x1, x2, params);

        % Compute errors
        diffp = p_comp - p_exact;
        pL2error(k) = sqrt(ds_x*ds_y*sum(diffp.^2));
        pLinferror(k) = max(max(abs(diffp)));
    end

    % Save Current Error results and update figure
    name = ['N = ' num2str(N)];
    
    figure(fig1); hold on;
    loglog(epvals, pL2error, 'o-', 'LineWidth', 2.5, 'MarkerSize', ...
       10,'DisplayName',name)
    hold off;

    figure(fig2); hold on;
    loglog(epvals, pLinferror,'x--','Linewidth',2.5, 'MarkerSize',10,'DisplayName',name)
    hold off;

    % Save Current Error results and update figure
    % with x-axis ep/ds
    name = ['N = ' num2str(N)];
    
    figure(fig3); hold on;
    semilogy(epvals/ds_x, pL2error,'o-', 'LineWidth', 2.5, 'MarkerSize', ...
       10,'DisplayName',name)
    hold off;

    figure(fig4); hold on;
    semilogy(epvals/ds_x,pLinferror,'x--','Linewidth',2.5, 'MarkerSize',10,'DisplayName',name)
    hold off;

end

figure(fig1)
box on
grid minor
set(gca,'YScale','log');
set(gca,'xscale','log');
legend('Location', 'Best')
fig_title = ['PressureError_L2_Permeable_StokesletsOnly','.fig'];
png_title = ['PressureError_L2_Permeable_StokesletsOnly','.png'];
savefig([currPath,'/',fig_title])
saveas(gcf, fullfile(currPath, png_title))

figure(fig2)
box on
grid minor
set(gca,'YScale','log');
set(gca,'xscale','log');
legend('Location', 'Best')
fig_title = ['PressureError_Linf_Permeable_StokesletsOnly','.fig'];
png_title = ['PressureError_Linf_Permeable_StokesletsOnly','.png'];
savefig([currPath,'/',fig_title])
saveas(gcf, fullfile(currPath, png_title))

figure(fig3)
box on
grid minor
set(gca,'YScale','log');
legend('Location', 'Best')
fig_title = ['PressureError_L2_Permeable_StokesletsOnly_epds','.fig'];
png_title = ['PressureError_L2_Permeable_StokesletsOnly_epds','.png'];
savefig([currPath,'/',fig_title])
saveas(gcf, fullfile(currPath, png_title))

figure(fig4)
box on
grid minor
set(gca,'YScale','log');
legend('Location', 'Best')
fig_title = ['PressureError_Linf_Permeable_StokesletsOnly_epds','.fig'];
png_title = ['PressureError_Linf_Permeable_StokesletsOnly_epds','.png'];
savefig([currPath,'/',fig_title])
saveas(gcf, fullfile(currPath, png_title))

% Save data for current run
clear fig1 fig2 fig3 fig4 % don't save figure handles
save(fullfile(currPath, 'RunData.mat'))