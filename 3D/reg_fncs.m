function [H1, H2, S1, S2, Q] = reg_fncs(ep, R, blob_num)

% Computes two components of blobs 
% for the Method of Regularized Stokeslets in 3D 
% Based on Cortez, Fauci, Medovikov, Physics of Fluids 2005 

% Developed by Shilpa Khatri and Ricardo Cortez, modified by Michaela
% Kubacki
% June 2024 

%ep: blob width (regularization parameter)
%R: distance between targe and source point + regularization 
%   (R = sqrt(|x-y|^2 + ep^2))

% Blob as given in Cortez, Fauci, Medovikov, Physics of Fluids 2005, Eqn 9 
A = (R.^2 + ep.^2) ./ (R.^3);
B = 1./ (R.^3);

% Another blob in 3D 
%R2 = R.^2; 
%R3 = R.^3; 
%d2 = ep.^2; 

%B = (2*R3.^2 +3*d2*R2.^2+35*d2^3)./(2*R3.^3) ;
%A = (R2 + d2).*B -2*d2*(15*d2^2+R2.^2)./(R3.*R2.^2) ;
r2 = R.^2 - ep^2;

switch blob_num
    case 1
        % blob = 15 d^4/(8*Pi*(r^2 + d^2)^(7/2))
        H1 = (ep^2 + R.^2) ./ (8 * pi * R.^3);
        H2 = 1 ./ (8*pi*R.^3);
        S1 = (3*ep^2 + 2*R.^2) ./ (8*pi*R.^5);
        S2 = -3*(5*ep^2 + 2*R.^2) ./ (8*pi*R.^7);
        Q = -105*ep^4 ./ (8*pi*R.^9);
   %case 2
        
end