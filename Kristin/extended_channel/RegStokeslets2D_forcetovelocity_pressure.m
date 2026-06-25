function [p] = RegStokeslets2D_forcetovelocity_pressure(y,f,x,ep,mu,blob_num,normal,wt)

% Computes velocities at set of points when given a set of points and
% forces at those points using the Method of Regularized Stokeslets 
% Based on Cortez, SIAM J. Sci Comput. 2001
% Constant flow is set to 0 here (Uo in Eqn 9)

% Developed by Shilpa Khatri and Ricardo Cortez July 2024. 
%
% This function requires access to function reg_funcs_withdoublet.m

% Modified December 2025 to incorporate a quad weight vector.
%
% Inputs
%       y = (y1,y2) source points
%       f = (f1,f2) forces at those source points (force density)
%       x = (x1,x2) target points 
%       u = (u1,u2) velocity evaluated at those target points 
%       mu = viscosity 
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

    norm1=normal(k,1); 
    norm2=normal(k,2); 
    
    normxy=norm1*XY1+norm2*XY2;

    fdotn = f1(k)*norm1 + f2(k)*norm2;
    
    p(:) = p(:) + fdotn*normxy.*S1*wt(k);

end

% rescaling
%u1 = u1/(mu); 
%u2 = u2/(mu); 

% repacking output 
%u = [u1 u2]; 


