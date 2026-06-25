function [f] = RegStokeslets2D_velocitytoforce_KK(y,x,u,ep,mu,blob_num,wt)

% Computes forces at set of points when given a set of points and
% velocities at those points using the Method of Regularized Stokeslets 
% Based on Cortez, SIAM J. Sci Comput. 2001

% Developed by Shilpa Khatri and Ricardo Cortez July 2024. Adapted by
% Michaela Kubacki December 2025.
% This function requires access to function reg_funcs_withdoublet.m
% Modified December 2025 to incorporate a quad weight vector.

% Inputs: 
%       y = (y1,y2) source points
%       f = (f1,f2) forces at those source points (force density) 
%       x = (x1,x2) target points 
%       u = (u1,u2) velocity at those target points 
%       mu is the viscosity 
%       ep is the width of the regularization 
%       blob_num = blob choice number
%       wt = quadrature weights for all source points 

% ACCURACTY NOTE: If inf norm of condition number greater than 10^r (r>0)
%                then losing r digits of accuracy in matrix solve. Warning
%                message will display if r >= 8.

N = size(y,1); %number of source points 

% unpacking velocity inputs
u1 = u(:,1);
u2 = u(:,2); 

% Assemble Stokeslet matrix
Mat = stokeslet_matrix(y,x,ep,mu,blob_num,wt);

% solving for force when given a velocity 
uu = [u1; u2];
ff = Mat\uu; 
f1 = ff(1:N); 
f2 = ff(N+1:end); 

% repacking output 
f = [f1, f2]; 

