function [p] = RegStokeslets2D_ftopressure(y, f, x, ep, mu, blob_num, wt)

% Computes the pressure on permeable from Stokeslet force f
% using the Method of Regularized Stokeslets Based on Cortez, SIAM
% J. Sci Comput. 2001 and Cortez, Fluids 2021

% Developed by Shilpa Khatri, Ricardo Cortez, Brittany Leathers, and
% Michaela Kubacki July 2024.
% Modified December 2025 to incorporate a quad weight vector.
% Modified June 2026 by Kurianski to compute pressure from stokeslets forces.

% This function requires access to function reg_funcs_withdoublet.m

% Inputs:
%       y = (y1,y2) source points 
%       f = (f1,f2) intermediate forces at those source points
%       x = (x1,x2) target points 
%       ep is the width of the regularization 
%       mu is the viscosity 
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

% initializing pressure
p = zeros(M,1);

% loop over source points    
for k = 1:N 

    % distance between target and source points 
    XY1 = x1(:) - y1(k); 
    XY2= x2(:) - y2(k);  
    
    % obtain source double solution pieces 
    Rsq = XY1.^2 + XY2.^2 + ep^2; 
    R2= sqrt( Rsq ); 

    [~, ~, S1, ~] = reg_fncs_withdoublet(ep,R2, blob_num); 

    % Calculate pressure 
    % p = sum_{k=1}^N (f_k * (x - y_k))S_1(r)
    fdotxy = f1(k)*XY1 + f2(k)*XY2;
    p(:) = p(:) + fdotxy.*S1*wt(k);

end