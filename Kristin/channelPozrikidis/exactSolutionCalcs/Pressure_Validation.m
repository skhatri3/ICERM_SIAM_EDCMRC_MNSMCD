% =========================================================================
% Validating the Pressure Solution
% =========================================================================
clear; clc; close all;

%% 1. Parameter Definitions
mu    = 1; % viscosity
Da    = 1; % Darcy coefficient
L     = 4; % Length of permeable region
c     = 2; % Length of extension region
G     = 1; % Pressure gradient
H     = 1; % Channel radius
alpha = 3*mu/(H^2*Da); % Related to permeability
gamma = -G/alpha * sech(alpha*L/2); % Constant from transmural pressure 
N_terms = 50; % Number of Fourier terms in the sum
W = L + 2*c; % Total domain width [0, W]
%% Validating the Fourier Coefficients for the top/bottom wall
x = linspace(0, W, 1000);

% Exact Piecewise Transmural Pressure p(x,0)=p(x,2H)= p_T(x)
p_exact = zeros(size(x));
for i = 1:length(x)
    xi = x(i);
    if xi >= 0 && xi < c
        p_exact(i) = gamma * sinh(alpha*(c - (L+2*c)/2)) - G*(xi - c);
    elseif xi >= c && xi <= (L + c)
        p_exact(i) = gamma * sinh(alpha*(xi - (L+2*c)/2));
    else % (L + c) < xi <= (L + 2*c)
        p_exact(i) = gamma * sinh(alpha*((L+c) - (L+2*c)/2)) - G*(xi - (L+c));
    end
end

% Fourier Coefficients for top/bottom problem
a0 = 0; % a0 = 0 by odd symmetry of p(x,0)
p_fourier = (a0 / 2) * ones(size(x));

for n = 1:N_terms
    lambda_n = n * pi / W; % eigenvalue from separation of variables    
    if mod(n, 2) == 0  % when n is EVEN (2, 4, 6...)
        % By odd symmetry of p(x,0) around W/2, integral over [0, W] is 0
        % This is because sinh term is odd and cos term is even.
        % so the integrand is odd.
        a_n = 0;        
    else % n is ODD (1, 3, 5...)
        % Now the integrand is even since cos term is odd
        % Integral over [0, c]
        An = -(gamma / lambda_n) * sinh(alpha*L / 2) * sin(lambda_n * c) ...
             + (G / (lambda_n^2)) * (1 - cos(lambda_n * c));
        % Integral over [c,L+2c] will be An by symmetry
      
        % Integral over [c, W/2]:
        Bn_half = (gamma / (alpha^2 + lambda_n^2)) * ...
          (-alpha * cosh(alpha*L / 2) * cos(lambda_n * c) ...
           + lambda_n * sinh(alpha*L / 2) * sin(lambda_n * c));
        % By symmetry, the other half [W/2,L+c] is the same.
               
        % Total integral over [0, W] is 2 * (Integral over [0, W/2])
        integral_total = 2 * (An + Bn_half);
        
        a_n = (2 / W) * integral_total;
    end    
    % Update Fourier Series sum
    p_fourier = p_fourier + a_n * cos(lambda_n * x);
end
%% Plotting Transmural Pressure Fourier Series
% Comparing the plot of continuous extension of p_T(x) to the Fourier series

figure('Color', 'w', 'Position', [100, 100, 900, 500]);
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
%% Visualizing the Full Pressure Solution
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

for k = 1:N_terms
    n = 2*k - 1;             % Odd mode: 1, 3, 5...
    lambda_n = n * pi / W;
    
    % Fourier coefficients a_n (for odd n)
    An = -(gamma / lambda_n) * sinh(alpha*L / 2) * sin(lambda_n * c) ...
         + (G / (lambda_n^2)) * (1 - cos(lambda_n * c));
     
    Bn_half = (gamma / (alpha^2 + lambda_n^2)) * ...
              (-alpha * cosh(alpha*L / 2) * cos(lambda_n * c) ...
               + lambda_n * sinh(alpha*L / 2) * sin(lambda_n * c));
           
    a_n = (4 / W) * (An + Bn_half);
    % Y_profile = (sinh(lambda_n*Y) - sinh(lambda_n*(Y-2H))) / sinh(2*lambda_n*H)
    % Using exp instead to avoid overflow errors
    Y_profile = exp(-lambda_n * (2*H - Y)) + exp(-lambda_n * Y);    
    p_tb = p_tb + a_n * Y_profile .* cos(lambda_n * X);
end

% Compute Left-Right Component: p_lr(x,y) (Overflow-Safe)
p_lr = zeros(size(X));

for k = 1:N_terms
    n = 2*k - 1;             % Odd mode: 1, 3, 5...
    mu_n = n * pi / (2*H);
    
    % (cosh(mu_n*(X-W)) - cosh(mu_n*X)) / sinh(mu_n*W)
    % Evaluates directly to [exp(-mu_n*X) - exp(-mu_n*(W-X))]
    % Using this to avoid overflow error.
    X_profile = exp(-mu_n * X) - exp(-mu_n * (W - X));
    
    % Coefficient without csch(mu_n*W) since it was absorbed above
    coeff = 1 / (n^2);
    
    p_lr = p_lr + coeff * X_profile .* sin(mu_n * Y);
end

p_lr = -(8 * G * H / (pi^2)) * p_lr;

% Total Pressure 
P_total = p_lr + p_tb;


% Visualization: Level Curves (Contours) + Colorbar
figure('Color', 'w', 'Position', [100, 100, 950, 600]);

% Filled color contour plot background
[~, h_fill] = contourf(X, Y, P_total, 35, 'LineStyle', 'none');
hold on;

% Explicit level curves overlay with contour labels
[C, h_lines] = contour(X, Y, P_total, 20, 'k-', 'LineWidth', 0.8);
clabel(C, h_lines, 'FontSize', 9, 'Color', 'k', 'LabelSpacing', 200);

% Domain markers for Region boundaries
xline(c, 'w--', 'LineWidth', 1.5, 'HandleVisibility', 'off');
xline(L+c, 'w--', 'LineWidth', 1.5, 'HandleVisibility', 'off');

% Formatting
colormap(turbo);           % Rich colormap for field intensity
cb = colorbar;
cb.Label.String = 'Pressure p(x,y)';
cb.Label.FontSize = 12;

xlabel('x', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('y', 'FontSize', 12, 'FontWeight', 'bold');
title('2D Pressure Field Solution p(x,y) with Contour Level Curves', 'FontSize', 14);

axis equal tight;
set(gca, 'FontSize', 11, 'Layer', 'top');
grid on;
