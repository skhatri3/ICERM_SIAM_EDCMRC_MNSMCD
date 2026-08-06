% =========================================================================
% Fourier Cosine Series Convergence for Piecewise Boundary Condition p(x,0)
% =========================================================================
clear; clc; close all;

%% 1. Parameter Definitions
L     = 4;          % Length of inner segment
c     = 2;           % Length of outer segments
G     = 1;         % Slope/constant parameter
alpha = 1;         % Hyperbolic decay/growth rate
gamma = -G/alpha * sech(alpha*L/2);
H     = 1;           % Domain height y
N_terms = 50;        % Number of Fourier terms in the sum

W = L + 2*c;         % Total domain width [0, W]

%% 2. Evaluation Grid
x = linspace(0, W, 1000);

%% 3. Exact Piecewise Function p(x,0)
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

%% 4. Fourier Cosine Series Reconstruction
a0 = 0; % a0 = 0 by odd symmetry of p(x,0)
p_fourier = (a0 / 2) * ones(size(x));

for n = 1:N_terms
    lambda_n = n * pi / W;
    
    if mod(n, 2) == 0  % n is EVEN (2, 4, 6...)
        % By odd symmetry of p(x,0) around W/2, integral over [0, W] is STRICTLY ZERO!
        a_n = 0;
        
    else               % n is ODD (1, 3, 5...)
        % Integral over Region 1 [0, c]:
        An = -(gamma / lambda_n) * sinh(alpha*L / 2) * sin(lambda_n * c) ...
             + (G / (lambda_n^2)) * (1 - cos(lambda_n * c));
         
        % Integral over Left Half of Region 2 [c, W/2]:
        Bn_half = (gamma / (alpha^2 + lambda_n^2)) * ...
          (-alpha * cosh(alpha*L / 2) * cos(lambda_n * c) ...
           + lambda_n * sinh(alpha*L / 2) * sin(lambda_n * c));
               
        % Total integral over [0, W] is 2 * (Integral over [0, W/2])
        integral_total = 2 * (An + Bn_half);
        
        a_n = (2 / W) * integral_total;
    end
    
    p_fourier = p_fourier + a_n * cos(lambda_n * x);
end
%% 5. Plotting Results
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


