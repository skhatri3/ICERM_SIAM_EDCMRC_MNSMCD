function [H1, H2, S1, S2, Q] = reg_fncs(ep, R, blob_num)

% Computes 5 components of blobs 
% for the Method of Regularized Stokeslets in 3D 
% Based on Cortez, Fauci, Medovikov, Physics of Fluids 2005 

% Developed by Shilpa Khatri and Ricardo Cortez, modified by Michaela
% Kubacki
% June 2024 

%ep: blob width (regularization parameter)
%R: distance between targe and source point + regularization 
%   (R = sqrt(|x-y|^2 + ep^2))

r2 = R.^2 - ep^2;

switch blob_num
    case 1
        % blob = 15 d^4/(8*Pi*(r^2 + d^2)^(7/2))
        H1 = (ep^2 + R.^2) ./ (8 * pi * R.^3);
        H2 = 1 ./ (8*pi*R.^3);
        S1 = (3*ep^2 + 2*R.^2) ./ (8*pi*R.^5);
        S2 = -3*(5*ep^2 + 2*R.^2) ./ (8*pi*R.^7);
        Q = -105*ep^4 ./ (8*pi*R.^9);
   case 2
        % blob = (15*d^4 (40*d^6 - 132*d^4*r^2 + 57*d^2*r^4 - 2*r^6))...
        %          /(2*(d^2 + r^2)^(13/2))
        H1 = (35*ep^8 - 25*ep^6*R.^2 + 3*ep^4*R.^4 + ep^2*R.^6 ...
            + 2*R.^8)./(2*R.^9);
        H2 = (35*ep^6 + 3*ep^2*R.^4 + 2*R.^6)./(2*R.^9);
        S1 = (315*ep^8 - 140*ep^6*R.^2 + 15*ep^4*R.^4 ...
            + 6*ep^2*R.^6 + 4*R.^8)./(2*R.^11);
        S2 = -3*(1155*ep^8 - 420*ep^6*R.^2 + 35*ep^4*R.^4 ...
            + 10*ep^2*R.^6 + 4*R.^8)./(2*R.^13);
        Q = 15*ep^4*(-3003*ep^6 + 3012*ep^4*R.^2 - 12*r2.^2*R.^2 ...
            + 26*R.^6 + 3*ep^2*(76*r2.*R.^2 - 273*R.^4))./(2*R.^15);
end