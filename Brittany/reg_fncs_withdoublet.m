function [H1, H2, S1, S2, Q] = reg_fncs_withdoublet(ep, R, blob_num)

% Computes two components of blobs (aka H1 and H2 in RC's notes) 
% for the Method of Regularized Stokeslets in 2D 
% Based on Cortez, SIAM J. Sci Comput. 2001

% Developed by Shilpa Khatri and Ricardo Cortez 
% July 2024 

%ep: blob width (regularization parameter)
%R: distance between targe and source point + regularization 
%   (R = sqrt(|x-y|^2 + ep^2))
%blob_num specifies blob type

switch blob_num
    case 1
        % phi = (2*d^4)/(pi*(r^2+d^2)^3) from 2021 paper
        % also called the "more common" blob
        % Min error roughly when ep = 0.95*ds
        H1 = (ep^2 - R.^2.*log(R))./(4*pi*R.^2);
        H2 = 1./(4*pi*R.^2); 
        S1 =(ep^2 + R.^2)./(2*pi*R.^4);
        S2 = -(2*ep^2 + R.^2)./(pi*R.^6);
        Q = - (12*ep^4)./(pi*R.^8);
    case 2
        % psi from Cortez Fluids 2021
        % psi = 2*ep^4*(r^4 - 10*ep^2*r^2 + 5*ep^4)/(Pi*(r^2 + ep^2)^5)  
        d2=ep^2;
        r2 = R.^2 - ep^2;
        H1 = (2/3*d2*(7*d2^2+r2.^2))./(8*pi*(r2+d2).^3) - log(r2+d2)/(8*pi);
        H2 = 2/3*(15*d2^2+10*d2*r2+3*r2.^2)./(8*pi*(r2+d2).^3);
        S1 = (10*d2^3 + 5*d2^2*r2 + 4*d2*r2.^2 + r2.^3)./(2*(r2+d2).^4)/pi;
        S2 = -(35*d2^3 + 7*d2^2*r2 + 5*d2*r2.^2 + r2.^3)./(  (r2+d2).^5)/pi;
        Q = (4*ep^4*(-80*ep^4 + 50*ep^2*R.^2) + 2*r2.*R.^2-5*R.^4)./(pi*R.^12);

end