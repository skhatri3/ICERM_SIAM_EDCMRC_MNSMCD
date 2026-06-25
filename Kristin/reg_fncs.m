function [H1, H2] = reg_fncs(ep, R, blob_num)

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
end