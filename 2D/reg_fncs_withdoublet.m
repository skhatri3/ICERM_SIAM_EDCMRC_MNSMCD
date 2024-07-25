function [H1, H2, S1, S2, Q] = reg_fncs_withdoublet(ep, R, blob_num)

% Computes two components of blobs (aka H1 and H2 in RC's notes) 
% for the Method of Regularized Stokeslets in 2D 
% Based on Cortez, SIAM J. Sci Comput. 2001

% Developed by Shilpa Khatri and Ricardo Cortez 
% July 2024 

%ep: blob width (regularization parameter)
%R: distance between targe and source point + regularization 
%   (R = sqrt(|x-y|^2 + ep^2))
%blob_num specifies blob type:
%   1 -- Blob as given in Cortez, SIAM J. Sci Comput. 2001, Eqns 11 
%   2 -- A more commonly used blob by RC 
%   3 -- Mystery blob from RC on 7/23

switch blob_num
    case 1
        %Blob as given in Cortez, SIAM J. Sci Comput. 2001, Eqns 11 
        H1 = -log(R + ep) + ep*(R + 2*ep)./(R + ep)./R; 
        H2 = (R + 2*ep)./(R + ep)./(R + ep)./R;
        r2=R.^2-ep^2;
        S1=-ep^3./(2*pi*r2.*R.^3);
        S2=ep^3*(2*ep^2+5*r2)./(2*pi*r2.^2.*R.^5);
        Q=-15*ep^3./(2*pi*R.^7);
    
    case 2
        %a more commonly used blob by RC 
        % Min error roughly when ep = 0.95*ds
        H1 = -log(R) +  ep^2./R./R;
        H2 = 1./R./R; 
    
    case 3
        % Mystery blob from RC on 7/23
        % Min error roughly when ep = 1.75*ds
        r2 = R.^2 - ep^2;
        H1 = (2/3*ep^2*(7*ep^4+r2.^2))./(r2+ep^2).^3 - log(r2+ep^2);
        H2 = 2/3*(15*ep^4+10*ep^2*r2+3*r2.^2)./(r2+ep^2).^3;
    case 4
        % psi from Cortez Fluids 2021
        % psi = 2*ep^4*(r^4 - 10*ep^2*r^2 + 5*ep^4)/(Pi*(r^2 + ep^2)^5)
        % From Mathematica:
        %   H1 = -((-11 d^6 + 9 d^4 r^2 + 7 d^2 r^4 + 3 r^6 + 3 (d^2 + r^2)^3 Log[d^2 + r^2])/(6 (d^2 + r^2)^3))
        %   H2 = (15 d^4 + 10 d^2 r^2 + 3 r^4)/(3 (d^2 + r^2)^3)
        %   S1 = (10 d^6 + 5 d^4 r^2 + 4 d^2 r^4 + r^6)/(2 \[Pi] (d^2 + r^2)^4)
        %   S2 = -((35 d^6 + 7 d^4 r^2 + 5 d^2 r^4 + r^6)/(\[Pi] (d^2 + r^2)^5))
        r2 = R.^2 - ep^2;
        H1 = -(-11*ep^6 + 9*ep^4*r2 + 7*d^2*r2.^2 + 3*r2.^3 + 3*R.^3.*log(R.^2))./(6*R.^6);
        H2 = (15*ep^4 + 10*ep^2*r2 + 3*r2.^2)./(3*R.^6);
        S1 = (10*ep^6 + 5*ep^4*r2 + 4*ep^2*r2.^2 + r2.^3)./(2*pi*R.^8);
        S2 = -(35*ep^6 + 7*ep^4*r2 + 5*ep^2*r2.^2 + r2.^3)./(pi*R^(10));
end