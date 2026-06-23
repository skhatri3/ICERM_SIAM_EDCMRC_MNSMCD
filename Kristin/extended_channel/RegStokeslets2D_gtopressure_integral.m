function [p, pjump] = RegStokeslets2D_gtopressure_integral(y,g,...
    x,ep,mu,blob_num,beta, normal, wt)

% Computes velocities on permeable membrane (u_beta) from intermediate
% force g using the Method of Regularized Stokeslets Based on Cortez, SIAM
% J. Sci Comput. 2001 and Cortez, Fluids 2021

% Developed by Shilpa Khatri, Ricardo Cortez, Brittany Leathers, and
% Michaela Kubacki July 2024.

% This function requires access to function reg_funcs_withdoublet.m

% Modified December 2025 to incorporate a quad weight vector.

% Inputs:
%       y = (y1,y2) source points 
%       g = (g1,g2) intermediate forces at those source points
%       x = (x1,x2) target points 
%       ep is the width of the regularization 
%       mu is the viscosity 
%       blob_num = blob choice number
%       beta = numerical permeability coefficient for all source points
%       beta(i) should be 0 if it is not a permeable section
%       normal = unit normal vectors for all source points
%       wt = quadrature weights for all source points 

% ************** NOTE ******************
% We are not using beta here. Should we be? We are using the
% Stokeslet pressure formula, with the f(s) = (g(s) dot n(s)) n(s). The
% pressure formulation we are using mirrors the formulation on p. 11 of the
% 2021 paper.

N = size(y,1); % number of source points 
M = size(x,1); % number of target points 

% unpacking the inputs 
y1 = y(:,1);
y2 = y(:,2);
g1 = g(:,1); 
g2 = g(:,2); 
x1 = x(:,1); 
x2 = x(:,2); 

% initializing the pressure and indicator variables
p = zeros(M,1);
chi = zeros(M,1); 

% For each target point, find closest source point 
close_idx = zeros(M,1); % for storing index of closest source point
for j = 1:M 
    d2_xy = (x1(j)-y1).^2 + (x2(j)-y2).^2; % distance squared for point xj and all source points
    [~, close_idx(j)] = min(d2_xy); % find the index of the closest source point 
    if abs(x2(j))<1 % if target point (x1,x2) is inside channel 
        chi(j) = -1; % indicator value is -1
    end % otherwise, chi(j) = 0
end

% g(s0) dot n(s0) at closest source points to target points (size of target points, Mx1)
gdotn_close = g1(close_idx).*normal(close_idx,1) + g2(close_idx).*normal(close_idx,2); 

% loop over source points    
for k = 1:N 
    % difference between target and source points 
    XY1 = x1(:) - y1(k); % Mx1
    XY2 = x2(:) - y2(k); % Mx1 
    
    % obtain stokeslet pressure solution pieces 
    Rsq = XY1.^2 + XY2.^2 + ep^2; 
    R2= sqrt( Rsq ); 

    [~, ~, S1, ~, ~] = reg_fncs_withdoublet(ep,R2, blob_num); 

    % Calculate pieces for stokeslet pressure formula
    n1=normal(k,1); 
    n2=normal(k,2); 
    
    nDotXY=n1*XY1+n2*XY2; % n dot (x-y)
    gdotn = g1(k)*n1 + g2(k)*n2;

    % Calculate stokeslet pressure
    
    p(:) = p(:) + (gdotn - gdotn_close).*nDotXY.*S1*wt(k);

end

% Add indicator piece of solution
pjump = gdotn_close.*chi;
p = p + pjump;





