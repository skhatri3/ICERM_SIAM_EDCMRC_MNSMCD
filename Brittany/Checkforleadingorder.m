clear all;
load C:\Users\britt\OneDrive\Documents\GitHub\ICERM_SIAM_EDCMRC_MNSMCD\Brittany\Error_eps_no_perm.mat
%%

[Nvals, epsvals]=ndgrid(Nvals, epsvals);
epsvals=reshape(epsvals, 27*5, 1);
Nvals=reshape(Nvals, 27*5,1);

E=reshape(eumaxnorm_away, 27*5,1);


Phi1 = [epsvals(:).^2, (1./Nvals(:)).^2./epsvals(:).^2];
Phi2 = [epsvals(:).^2, (1./Nvals(:)).^2./epsvals(:)];
Phi3 = [epsvals(:).^2, (1./Nvals(:)).^2];
Phi4 = [epsvals(:).^2, (1./Nvals(:)).^2.*abs(log(epsvals(:)))];

coef1 = Phi1\E(:) 
res1 = norm(Phi1*coef1-E(:))/norm(E(:))
coef2 = Phi2\E(:)  
res2 = norm(Phi2*coef2-E(:))/norm(E(:))
coef3 = Phi3\E(:)  
res3 = norm(Phi3*coef3-E(:))/norm(E(:))
coef4 = Phi4\E(:) 
res4 = norm(Phi4*coef4-E(:))/norm(E(:))