function [u_beta] = RegStokeslets2D_gtovelocity(y,g,...
    x,ep,mu,blob_num, beta, normal, wt)

% Computes velocities on permeable membrane (u_beta) from intermediate
% force g using the Method of Regularized Stokeslets Based on Cortez, SIAM
% J. Sci Comput. 2001 and Cortez, Fluids 2021

% Developed by Shilpa Khatri, Ricardo Cortez, Brittany Leathers, and
% Michaela Kubacki July 2024.

% This function requires access to function reg_funcs_withdoublet.m

% Modified December 2025 to incorporate a quad weight vector and enforce
% consistent scaling.

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

N = size(y,1); % number of source points 
M = size(x,1); % number of target points 

% unpacking the inputs 
y1 = y(:,1);
y2 = y(:,2);
g1 = g(:,1); 
g2 = g(:,2); 
x1 = x(:,1); 
x2 = x(:,2); 

% initializing the velocity 
u1 = zeros(M,1);
u2 = zeros(M,1); 

% loop over source points    
for k = 1:N 

    % distance between target and source points 
    XY1 = x1(:) - y1(k); 
    XY2= x2(:) - y2(k);  
    
    % obtain source double solution pieces 
    Rsq = XY1.^2 + XY2.^2 + ep^2; 
    R2= sqrt( Rsq ); 

    [~, ~, S1, S2, ~] = reg_fncs_withdoublet(ep,R2, blob_num); 

    % Calculate u 
    norm1=normal(k,1); 
    norm2=normal(k,2); 
    
    normxy=norm1*XY1+norm2*XY2;

    u1(:)=u1(:)-(beta(k)*norm1*(S1.*norm1+S2.*normxy.*XY1)*g1(k)+...
        -beta(k)*norm2*(S2.*normxy.*XY1)*g2(k))*wt(k);

    u2(:)=u2(:)-(beta(k)*norm2*(S1.*norm2+S2.*normxy.*XY2)*g2(k)+...
        -beta(k)*norm1*(S2.*normxy.*XY2)*g1(k))*wt(k);    

end

% rescaling
u1 = u1/(8*pi*mu); 
u2 = u2/(8*pi*mu); 

% repacking output 
u_beta = [u1 u2]; 



