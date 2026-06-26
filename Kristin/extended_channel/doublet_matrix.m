function [Mat] = doublet_matrix(y,x,ep,mu,blob_num,wt,beta,normals)

% Computes doublet matrix for Method of Regularized Stokeslets incoporating
% permeability 
% Used to computes forces at set of points when given a set of points and
% velocities at those points using the Method of Regularized Stokeslets 
% Based on Cortez, SIAM J. Sci Comput. 2001 and Cortez, Fluids 2021 

% Developed by Shilpa Khatri and Ricardo Cortez June 2026

%y = (y1,y2) source points Nx2
%x = (x1,x2) target points Mx2 
%mu is the viscosity 
%ep is the width of the regularization for source doublet
%wt is the weights for the quadrature Nx1
%beta gives the permeability coefficient for all source points Nx1
%normals is Nx2, giving unit normals for source points

N = size(y,1); %number of source points 
M = size(x,1); %number of target points 

%unpacking the inputs 
y1 = y(:,1);
y2 = y(:,2);
x1 = x(:,1); 
x2 = x(:,2); 
n1 = normals(:,1); 
n2 = normals(:,2); 