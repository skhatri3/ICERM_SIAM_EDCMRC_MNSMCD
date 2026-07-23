function [p,pst,psd] = RegStokeslets2D_gtopressure_diff_ep(y,g,...
    x,ep1,ep2,mu,blob_num, beta, normal, wt)

% Computes the pressure on permeable membrane (u_beta) from intermediate
% force g using the Method of Regularized Stokeslets Based on Cortez, SIAM
% J. Sci Comput. 2001 and Cortez, Fluids 2021

% Developed by Shilpa Khatri, Ricardo Cortez, Brittany Leathers, and
% Michaela Kubacki July 2024.

% This function requires access to function reg_funcs_withdoublet.m

% Modified December 2025 to incorporate a quad weight vector.

% ************** NOTE ******************
% We are not using beta here. Should we be? We are using the
% Stokeslet pressure formula, with the f(s) = (g(s) dot n(s)) n(s).

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
pst = zeros(M,1);
psd = zeros(M,1);

% loop over source points    
for k = 1:N 

    % distance between target and source points 
    XY1 = x1(:) - y1(k); 
    XY2= x2(:) - y2(k);  
    
    % obtain source double solution pieces

    % Stokeslet Piece
    Rsq = XY1.^2 + XY2.^2 + ep1^2; 
    R1= sqrt( Rsq ); 

    [~, ~, S1, ~, ~] = reg_fncs_withdoublet(ep1,R1, blob_num); 
    
    % Source Doublet Piece
    Rsq = XY1.^2 + XY2.^2 + ep2^2; 
    R2= sqrt( Rsq ); 

    [~, ~, ~, ~, Q] = reg_fncs_withdoublet(ep2,R2, blob_num);
    
    % Calculate p 
    norm1=normal(k,1); 
    norm2=normal(k,2); 
    
    ndotxy=norm1*XY1+norm2*XY2;
    gdotn = g1(k)*norm1 + g2(k)*norm2;
    
    % Stokeslet pressure
    pst = pst(:) + beta(k)*gdotn*ndotxy.*S1*wt(k);

    % Source Doublet Pressure
    psd = psd(:) - (beta(k)*gdotn)*ndotxy.*Q*wt(k);
     
end

p = pst + psd;



