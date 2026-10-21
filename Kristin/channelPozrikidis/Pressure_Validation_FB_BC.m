% =========================================================================
% Fourier Cosine Series Convergence for Piecewise Boundary Condition
% p(x,0)=p(x,2H)
% =========================================================================
clear; clc; close all;

%% Parameters
L = 3;    % Length of permeable region
c = 2;    % Length of extension region
G = 4;    % Pressure gradient
Da = 0.4; % Darcy number
H = 1;    % Channel radius
xm = (L+2*c)/2; % Channel midpoint

% Use Newton's method to compute lambda eigenvalue
eta = sqrt(2*Da);
fun = @(L) (Da-H/2)*L*cot(L*H)*cos(L*H) + cos(L*H)/2 - H*L*sin(L*H)/2;
eta = fzero(fun, eta);

bL = G*(c + tanh(eta*L/2)/eta);
bR = G*(L+2*c) - bL;

N_vec = 50;%[50, 100, 200, 400]; % Number of Fourier terms in the sum
W = L + 2*c;         % Total domain width [0, W]

x = linspace(0, W, 1000);

%% Construct exact piecewise function for p(x,0)=p(x,2H)
% For use on top/bottom boundary
p_exact = zeros(size(x));

for i = 1:length(x)
    x_i = x(i);
    if x_i >= 0 && x_i < c
        p_exact(i) = G*(c - x_i + tanh(eta*L/2)/eta);
    elseif x_i >= c && x_i <= (L + c)
        p_exact(i) = -G*sinh(eta*(x_i - xm))/(eta*cosh(eta*L/2));
    else % (L + c) < xi <= (L + 2*c)
        p_exact(i) = G*(L + c - x_i - tanh(eta*L/2)/eta);
    end
end

% Plot transmural pressure
figure
colors = get(gca, 'ColorOrder'); % get matlab default plotting colors
plot([c,c],[-15,15],'--k',[L+c,L+c],[-15,15], '--k','LineWidth',1.3)
hold on
plot([0,L+2*c],[0,0], '-k', [c+L/2, c+L/2], [-0.5,0.5], '-k', 'LineWidth',1.3)
plot(x,p_exact,'linewidth',2,'Color',colors(1,:))
% format
grid on
box on
text(c-0.3,-0.6,'$c$', 'Interpreter','latex','FontSize', 20)
text(L+c+0.1,-1,'$L+c$', 'Interpreter','latex','FontSize', 20)
text(c+L/2,-1.5,'$x_m$', 'Interpreter','latex','FontSize', 20)
xlabel('$x$', 'Interpreter', 'latex')
ylabel('$p_T(x)$', 'Interpreter', 'latex')
ax = gca; ax.FontSize = 20;
xlim([0, L+2*c])

%% Compute Fourier cosine series
maxErrAbs = zeros(size(N_vec));
L2Err = zeros(size(N_vec));
id = 1;
for N_terms = N_vec
    p_fourier = zeros(size(x));
    for n = 1:2:N_terms % loop over odd terms
        lambda_n = n*pi/W;

        % Fourier coefficients
        an = (eta^2*(1 - cos(c*lambda_n)) + lambda_n^2 + eta*lambda_n*sin(c*lambda_n)*tanh(L*eta/2))/(lambda_n^2*(eta^2 + lambda_n^2));

        % Update Fourier Series sum
        p_fourier = p_fourier + an*cos(lambda_n*x);
    end
    p_fourier = (4*G/W)*p_fourier; % multiply by constant from Fourier computation


    % Error computation
    errAbs = abs(p_fourier - p_exact);
    maxErrAbs(id) = max(errAbs);
    L2Err(id) = norm(p_fourier - p_exact);
    id = id+1;

    % Plot results
    % error plots
    figure(10)
    hold on
    plot(x, errAbs, 'LineWidth', 1.3)
end

%formatting
title('Fourier cosine series error for pressure boundary condition')
box on
xlabel('$x$', 'Interpreter', 'latex')
ylabel('Absolute error')
if length(N_vec) == 4
    legend(['N=', num2str(N_vec(1)) ', max abs error=', num2str(maxErrAbs(1)), ', L2 error=', num2str(L2Err(1))],...
        ['N=', num2str(N_vec(2)) ', max abs error=', num2str(maxErrAbs(2)), ', L2 error=', num2str(L2Err(2))],...
        ['N=', num2str(N_vec(3)) ', max abs error=', num2str(maxErrAbs(3)), ', L2 error=', num2str(L2Err(3))],...
        ['N=', num2str(N_vec(4)) ', max abs error=', num2str(maxErrAbs(4)), ', L2 error=', num2str(L2Err(4))],...
        'Location', 'best')
end
ax = gca; ax.FontSize = 14;

% plot of exact and series
figure
plot(x, p_exact, 'r-', 'LineWidth', 2.5, 'DisplayName', 'Exact Piecewise p(x,0)');
hold on;
plot(x, p_fourier, 'b--', 'LineWidth', 1.5, 'DisplayName', sprintf('Fourier Cosine Series (N = %d)', N_terms));

% Formatting
grid on;
xline(c, 'k:', 'LineWidth', 1.2, 'HandleVisibility', 'off');
xline(L+c, 'k:', 'LineWidth', 1.2, 'HandleVisibility', 'off');
xlabel('x', 'FontSize', 12);
ylabel('p(x,0)', 'FontSize', 12);
title('Convergence of Fourier Cosine Series to Piecewise Boundary Condition', 'FontSize', 14);
legend('Location', 'northwest', 'FontSize', 11);
set(gca, 'FontSize', 11);


%% Linearization approach
% Treat hyperbolic sine as a perturbation to linear
%xi = x(x>c & x<L+c);
p_linear = bL*(1-2*x/(L+2*c));
ppert_exact = p_exact - p_linear;

figure
plot(x, p_exact, 'LineWidth', 2)
hold on
plot(x, p_linear, '--', 'LineWidth', 2)
plot(x, ppert_exact, 'k', 'LineWidth', 1.3)

legend('Hyp. sine pressure', 'Constant press. grad.', 'Difference')
grid on
xlim([0, L+2*c])

%% Check Fourier cosine series for ppert_exact
% Compute Fourier cosine series
maxErrAbs = zeros(size(N_vec));
L2Err = zeros(size(N_vec));
id = 1;
for N_terms = N_vec
    ppert_fourier = zeros(size(x));
    for n = 1:2:N_terms % loop over odd terms
        lambda_n = n*pi/W;

        % Fourier coefficients
        an = (eta^2*(1 - cos(c*lambda_n)) + lambda_n^2 + eta*lambda_n*sin(c*lambda_n)*tanh(L*eta/2))/(lambda_n^2*(eta^2 + lambda_n^2));

        % Update Fourier Series sum
        ppert_fourier = ppert_fourier + an*cos(lambda_n*x);
    end
    ppert_fourier = (4*G/W)*ppert_fourier; % multiply by constant from Fourier computation


    % Error computation
    errAbs = abs(ppert_fourier - ppert_exact);
    maxErrAbs(id) = max(errAbs);
    L2Err(id) = norm(ppert_fourier - ppert_exact);
    id = id+1;

    % Plot results
    % error plots
    figure(11)
    hold on
    plot(x, errAbs, 'LineWidth', 1.3)
end

%formatting
title('Fourier cosine series error for pressure boundary condition')
box on
xlabel('$x$', 'Interpreter', 'latex')
ylabel('Absolute error')
if length(N_vec) == 4
    legend(['N=', num2str(N_vec(1)) ', max abs error=', num2str(maxErrAbs(1)), ', L2 error=', num2str(L2Err(1))],...
        ['N=', num2str(N_vec(2)) ', max abs error=', num2str(maxErrAbs(2)), ', L2 error=', num2str(L2Err(2))],...
        ['N=', num2str(N_vec(3)) ', max abs error=', num2str(maxErrAbs(3)), ', L2 error=', num2str(L2Err(3))],...
        ['N=', num2str(N_vec(4)) ', max abs error=', num2str(maxErrAbs(4)), ', L2 error=', num2str(L2Err(4))],...
        'Location', 'best')
end
ax = gca; ax.FontSize = 14;

% plot of exact and series
figure
plot(x, ppert_exact, 'r-', 'LineWidth', 2.5, 'DisplayName', 'Exact Piecewise p(x,0)');
hold on;
plot(x, ppert_fourier, 'b--', 'LineWidth', 1.5, 'DisplayName', sprintf('Fourier Cosine Series (N = %d)', N_terms));

% Formatting
grid on;
xline(c, 'k:', 'LineWidth', 1.2, 'HandleVisibility', 'off');
xline(L+c, 'k:', 'LineWidth', 1.2, 'HandleVisibility', 'off');
xlabel('x', 'FontSize', 12);
ylabel('p(x,0)', 'FontSize', 12);
title('Convergence of Fourier Cosine Series to Piecewise Boundary Condition', 'FontSize', 14);
legend('Location', 'northwest', 'FontSize', 11);
set(gca, 'FontSize', 11);








%% =========================================================================
% 2D Pressure Field Reconstruction via Fourier Superposition
% Solves Laplace's Equation: p(x,y) = p_lr(x,y) + p_tb(x,y)
% =========================================================================
%
% p_lr(x,y) is the solution to the Left/Right BCs 
% 
%       p_x(0,y) = p_x(L+2c,y) = -G, p(x,0)=p(x,2H) = 0
%
% p_tb(x,y) is the solution to the Top/Bottom BCs
%
%       p_x(0,y) = p_x(L+2c,y) = 0, p(x,0) = p(x,2H) = p_T(x)

% Grid Setup
Nx = 200; 
Ny = 200;
x = linspace(0, W, Nx);
y = linspace(0, 2*H, Ny);
[X, Y] = meshgrid(x, y);

% Compute Top-Bottom Component: p_tb(x,y) (Overflow-Safe)
p_tb = zeros(size(X));

for n = 1:2:N_terms % loop over odd terms
    lambda_n = n * pi / W; % Note: n = 2k-1
    
    % Fourier coefficients a_n (for odd n)
    a_n = (eta^2*(1 - cos(c*lambda_n)) + lambda_n^2 + eta*lambda_n*sin(c*lambda_n)*tanh(L*eta/2))/(lambda_n^2*(eta^2 + lambda_n^2));
    
    Y_profile = exp(-lambda_n*(2*H - Y)) + exp(-lambda_n * Y);
    
    p_tb = p_tb + a_n*Y_profile.*cos(lambda_n * X);
end
p_tb = (4*G/W)*p_tb; % multiply by constant in Fourier computation

% Compute Left-Right Component: p_lr(x,y) (Overflow-Safe)
p_lr = zeros(size(X));

for n = 1:2:N_terms % loop over odd terms
    mu_n = n * pi / (2*H);
    
    % Numerically stable horizontal ratio:
    % (cosh(mu_n*(X-W)) - cosh(mu_n*X)) / sinh(mu_n*W)
    % Evaluates directly to [exp(-mu_n*X) - exp(-mu_n*(W-X))]
    X_profile = exp(-mu_n * X) - exp(-mu_n * (W - X));
    
    % Coefficient without csch(mu_n*W) since it was absorbed above
    coeff = 1 / (n^2);
    
    p_lr = p_lr + coeff * X_profile .* sin(mu_n * Y);
end

p_lr = -(8*G*H/(pi^2))*p_lr;

% Total Pressure Field
P_total = p_lr + p_tb;


% Visualization: Level Curves (Contours) + Colorbar
figure('Color', 'w', 'Position', [100, 100, 950, 600]);

% Filled color contour plot background
[~, h_fill] = contourf(X, Y, P_total, 35, 'LineStyle', 'none');
hold on;
% level curves overlay with contour labels
[C, h_lines] = contour(X, Y, P_total, 20, 'k-', 'LineWidth', 0.8);
clabel(C, h_lines, 'FontSize', 14, 'Color', 'k', 'LabelSpacing', 200);
% Domain markers for region boundaries
xline(c, 'w--', 'LineWidth', 1.5, 'HandleVisibility', 'off');
xline(L+c, 'w--', 'LineWidth', 1.5, 'HandleVisibility', 'off');

% formatting
colormap(parula);           % Rich colormap for field intensity
cb = colorbar;
cb.Label.String = 'Pressure $p(x,y)$';
cb.Label.Interpreter = 'latex';
cb.Label.FontSize = 20;
xlabel('$x$', 'FontSize', 20, 'Interpreter','latex');
ylabel('$y$', 'FontSize', 20, 'Interpreter','latex');
%title('2D Pressure Field Solution p(x,y) with Contour Level Curves', 'FontSize', 14);
axis equal tight;
set(gca, 'FontSize', 20, 'Layer', 'top');
grid on;