clear all;
load C:\Users\britt\OneDrive\Documents\GitHub\ICERM_SIAM_EDCMRC_MNSMCD\Brittany\Error_eps_no_perm.mat
%%

[Nvals, epsvals]=ndgrid(Nvals, epsvals);
epsvals=reshape(epsvals, 27*5, 1);
Nvals=reshape(Nvals, 27*5,1);

Phi = [epsvals(:).^2, ...
    (1./Nvals(:)).^2 ./ epsvals(:).^2, ...
    (1./Nvals(:)).^2 ./ epsvals(:), ...
    (1./Nvals(:)).^2 .* abs(log(epsvals(:)))];
%%
coef = Phi \ reshape(eumaxnorm_near, 5*27, 1)