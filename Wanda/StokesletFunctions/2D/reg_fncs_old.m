function [A, B] = reg_fncs(ep,R)

% Computes two components of blobs (aka H1 and H2 in RC's notes) 
% for the Method of Regularized Stokeslets in 2D 
% Based on Cortez, SIAM J. Sci Comput. 2001

% Developed by Shilpa Khatri and Ricardo Cortez 
% July 2024 

%ep: blob width (regularization parameter)
%R: distance between targe and source point + regularization 
%   (R = sqrt(|x-y|^2 + ep^2))

%Blob as given in Cortez, SIAM J. Sci Comput. 2001, Eqns 11 
%A = -log(R + ep) + ep*(R + 2*ep)./(R + ep)./R; 
%B = (R + 2*ep)./(R + ep)./(R + ep)./R;

%a more commonly used blob by RC
A = -log(R) +  ep^2./R./R;
B = 1./R./R;

% r2 = R.^2-ep^2;
% A = (2/3*ep^2*(7*ep^4+r2.^2))./(r2+ep^2).^3 - log(r2+ep^2); % H1
% B = 2/3*(15*ep^4+10*ep^2*r2+3*r2.^2)./(r2+ep^2).^3; % H2
% 
% A = A/2;
% B = B/2;

   
