function [p] = RegStokeslets2D_forcetopressure_gauss(y,f,x,ep,mu,beta,blob_num,normal,tangent,wt)

% Computes pressure at set of points when given a set of points and %
% forces at those points using the pressure formula from the Method of
% Regularized Stokeslets after implementing Gauss's Law. 
% 
% Based on pressure formulation from Cortez et. al Fluids 2021, Example 3.

% This adaptation is for the infinite channel problem. The indicator
% function, chi, is specific to the geometry of this problem.

% Developed by Michaela Kubacki June 2026. 
%
% This function requires access to function reg_funcs_withdoublet.m
%
% Inputs
%       y = (y1,y2) source points
%       f = (f1,f2) forces at those source points (force density)
%       x = (x1,x2) target points 
%       u = (u1,u2) velocity evaluated at those target points 
%       mu = viscosity 
%       normal = normal vectors for the entire boundary
%       tau = tangent vectors for the entire boundary
%       ep is the width of the regularization 
%       blob_num = blob choice number
%       wt = quadrature weights for all source points 

N = size(y,1); % number of source points 
M = size(x,1); % number of target points 

% unpacking the inputs 
y1 = y(:,1);
y2 = y(:,2);
f1 = f(:,1);
f2 = f(:,2); 
x1 = x(:,1); 
x2 = x(:,2);
wt = wt(:); % making sure wt is a column vector

% initializing the pressure and indicator variables
p = zeros(M,1);
chi = zeros(M,1); 

% For each target point, find closest source point 
close_idx = zeros(M,1); % for storing index of closest source point

for j = 1:M % looping over target points
    d2_xy = (x1(j)-y1).^2 + (x2(j)-y2).^2; % d^2 for point xj and all source points (Nx1)
    [~, close_idx(j)] = min(d2_xy); % find the index of the closest source point 
    if abs(x2(j))<1 % if target point (x1,x2) is inside channel 
        chi(j) = -1; % indicator value is -1
    end % otherwise, chi(j) = 0
end

% Evaluate f dot n at the closest source point to each target point
% Size is Mx1
fdotn_close = f1(close_idx).*normal(close_idx,1)+f2(close_idx).*normal(close_idx,2);

% initializing the velocity 
%u1 = zeros(M,1);
%u2 = zeros(M,1);
p = zeros(M,1);

% loop over source points    
for k = 1:N 

    % distance between target and source points 
    XY1 = x1(:) - y1(k); 
    XY2 = x2(:) - y2(k);  
    R2 = XY1.^2 + XY2.^2 + ep^2; 
    R = sqrt( R2 ); 

    % computing the velocity 
    [~, ~, S1, ~] = reg_fncs_withdoublet(ep,R,blob_num);

    n1=normal(k,1); 
    n2=normal(k,2);

    tau1 = tangent(k,1);
    tau2 = tangent(k,2);

    nDotXY = n1*XY1 + n2*XY2; % n dot (x-y)
    tauDotXY = tau1*XY1 + tau2*XY2;

    fdotn = f1(k)*n1 + f2(k)*n2;

    fdottau = f1(k)*tau1 + f2(k)*tau2;

    % Calculate Stokeslet pressure
    
    p(:) = p(:) + (fdotn - fdotn_close).*nDotXY.*S1*wt(k) + fdottau.*tauDotXY.*S1*wt(k);

end

% Add indicator piece of solution
p = p + fdotn_close.*chi;
