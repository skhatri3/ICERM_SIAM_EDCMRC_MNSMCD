function p_fourier = exactPressureFun_Fourier(X, Y, params)

% INPUTS:
% x1gg -- coordinate grid for x1
% x2gg -- coordinate grid for x2
% params -- parameters for geometry and set-up

% Unpack parameters
G = params.G;
eta = params.eta;
L = params.L;
c = params.c;
H = params.H;

% Computes exact pressure in channel using Fourier series construction
bL = G*(c + tanh(eta*L/2)/eta); % note: bR = G*(L+2*c) - bL;

N_terms = 100;%[50, 100, 200, 400]; % Number of Fourier terms in the sum
W = L + 2*c;         % Total domain width [0, W]

% Compute Top-Bottom Component: p_tb(x,y) (Overflow-Safe)
p_tb = zeros(size(X));
for n = 1:2:N_terms % loop over odd terms
    lambda_n = n*pi / W; % Note: n = 2k-1
    
    % Fourier coefficients a_n (for odd n)
    a_n = (eta^2*(1 - cos(c*lambda_n)) + lambda_n^2 + eta*lambda_n*sin(c*lambda_n)*tanh(L*eta/2))/(lambda_n^2*(eta^2 + lambda_n^2));
    
    Y_profile = (exp(-lambda_n*Y) - exp(-lambda_n*(4*H-Y)) + exp(-lambda_n*(2*H-Y)) ...
        - exp(-lambda_n*(2*H+Y)))/(1 - exp(-4*H*lambda_n));
    
    p_tb = p_tb + a_n*Y_profile.*cos(lambda_n * X);
end
p_tb = (4*G/W)*p_tb; % multiply by constant in Fourier computation

% Compute Left-Right Component: p_lr(x,y) (Overflow-Safe)
p_lr = zeros(size(X));
for n = 1:2:N_terms % loop over odd terms
    mu_n = n*pi / (2*H);
    
    % Numerically stable horizontal ratio:
    % (cosh(mu_n*(X-W)) - cosh(mu_n*X)) / sinh(mu_n*W)
    X_profile = (exp(-mu_n*X) - exp(-mu_n*(W - X)))/(1 + exp(-mu_n*W));
    
    % Coefficient without csch(mu_n*W) since it was absorbed above
    p_lr = p_lr + (1/(n^2))*X_profile.* sin(mu_n*Y);
end

p_lr = -(8*G*H/(pi^2))*p_lr;

% Total Pressure Field
p_fourier = (p_lr + p_tb);


end