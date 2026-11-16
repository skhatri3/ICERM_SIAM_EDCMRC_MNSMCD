%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Developed for ICERM
% Exact solution to 3D Stokeslets
% Cortez, Fauci, Medovikov 2005 Section IV: Example 1
% Date: June 25, 2026
% Author: Wanda Strychalski
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all
clear all
clc
%% Parameters to set
point_num  =3;
factor = 1;
blob_num = 1;
mu = 1;       % fluid viscosity
%alpha = 1e-1; % permeability parameter
r = 1;

[xf,yf,zf] = TdesignW(point_num);
pts = [xf(:),yf(:),zf(:)];

X = pts(:,1);
Y = pts(:,2);
Z = pts(:,3);
ptN = length(X);

dA0 = 4*pi/ptN;       % try using a uniform weight
pointArea = dA0*ones(ptN,1);% try using a uniform weight
ep = 2*factor*sqrt(dA0); % epsilon or delta for MRS

U = 1;
Uexact = zeros(ptN,3);
Uexact(:,3) = U;
f = zeros(ptN,3);
f(:,3) = 3*mu/(2*r)*U;

z = [X,Y,Z];
f_weighted = f * dA0;
u = RegStokeslets3D_forcetovelocity(z,f_weighted,z,ep,mu,blob_num);
error = abs(u-Uexact);
for j=1:ptN
error_inf(j) = norm(error(j,:), Inf);
end

error_norm = norm(error_inf,'inf');

figure(1)
scatter3(X,Y,Z,'filled');
hold on
quiver3(X,Y,Z,u(:,1),u(:,2),u(:,3),'-k');
axis equal
%length = norm()
hold off
