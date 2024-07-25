function [H1, H2, S1, S2] = reg_fncs_withdoublet(ep, R, blob_num)

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
        %   H1 = (8 d^6 - 2 d^4 R^2 + d^2 R^4 - 3 R^6 Log[R])/(3 R^6)
        %   H2 = (8 d^4 + 4 d^2 R^2 + 3 R^4)/(3 R^6)
        %   S1 = (8 d^6 + d^2 R^4 + R^6)/(2 \[Pi] R^8)
        %   S2 = -((32 d^6 + 2 d^2 R^4 + R^6)/(\[Pi] R^10))
        H1 = (8*ep^6 - 2*ep^4*R.^2 + ep^2*R.^4 - 3*R.^6*log(R))./(3*R.^6);
        H2 = (8*ep^4 + 4*ep^2*R.^2 + 3*R.^4)./(3*R.^6);
        S1 = (8*ep^6 + ep^2*R.^4 + R.^6)./(2*pi*R.^8);
        S2 = -((32*ep^6 + 2*ep^2*R.^4 + R.^6)./(pi*R.^10));
end