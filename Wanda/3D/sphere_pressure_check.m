close all
clear all
clc
addpath('../StokesletFunctions/3D/');

%% Parameters to set
point_num  = 3;
factor = 1;
blob_num = 1;
mu = 1;       % fluid viscosity

r = 1;

[xf,yf,zf] = TdesignW(point_num);
pts = r*[xf(:),yf(:),zf(:)];

X = pts(:,1);
Y = pts(:,2);
Z = pts(:,3);
ptN = length(X);
N = ptN;

% Create a perfect triangulation topology based on the uniform sphere
surface_faces = boundary(X, Y, Z);
% Wrap in a triangulation object to easily compute face/vertex normals and areas
TR = triangulation(surface_faces, X, Y, Z);

dA0 = 4*pi*r^2/ptN;       % try using a uniform weight
pointArea = dA0*ones(ptN,1);% try using a uniform weight

ep = 2*factor*sqrt(dA0); % epsilon or delta for MRS

% % For the permeable code
% if blob_num==1 
%     beta = 4/3*ep*alpha*ones(ptN,1);
% elseif blob_num == 3
%     % ep = ep/2.1545; % to match the height of blob 1
%     beta = 4/9*ep*alpha*ones(ptN,1);
% else
%     beta = 0;
% end

% Convert to spherical coordinate to exactly 
% compute the normal vectors and curvature
rho = sqrt(X.^2 + Y.^2 + Z.^2);
theta = atan2(Y,X);
phi = acos(Z./rho);

exact_normals = [cos(theta).*sin(phi),sin(theta).*sin(phi),cos(phi)];
quiver3(X,Y,Z,exact_normals(:,1),exact_normals(:,2),exact_normals(:,3));
axis equal

R = 1/N*sum(sqrt(X.^2 + Y.^2+Z.^2));
if(abs(R)<= 1e-2)
    R = 1e-2;
    kappa_vec = 1e2*ones(N,1);
else
    kappa_vec = 1./(sqrt(X.^2 + Y.^2+Z.^2)); % curvature at each point as 1/R.
end
kappa_vec = 2.0*kappa_vec; % mean curvature

forces = (kappa_vec).*exact_normals; % checking with the exact normals.

dx = 3/50;
x = [-1.5:dx:1.5];
nx = length(x);
y = [-1.5:dx:1.5];
ny = length(y);
z = [-1.5:dx:1.5];
nz = length(z);
[X1,Y1,Z1] = meshgrid(x,y,z);
Xq =  [X1(:),Y1(:),Z1(:)];
% Reshape into an N x 3 matrix of [X, Y, Z] points

xT = [X,Y,Z];

forces = -forces*dA0;
 %[p] = RegStokeslets3D_forcetopressure(xT,forces,Xq,ep,blob_num);
  [p] = RegStokeslets3D_forcetopressure_fix(xT, kappa_vec, Xq, ep, blob_num,...
    exact_normals, pointArea, surface_faces);

 % 1. Reshape p back into the 3D grid dimensions [ny, nx, nz]
P3D = reshape(p, ny, nx, nz);
pmin = min(P3D(:));
deltaP = max(P3D(:)) - pmin
%P3D = P3D-pmin;

% 2. Select the index along the y-axis for the XZ-plane slice
% Find the index corresponding to y closest to 0
[~, y_idx] = min(abs(y - 0)); 

% 3. Extract the 2D slice across all x and z for that y-index
% P_slice will have dimensions [nx, nz]
P_slice = squeeze(P3D(y_idx, :, :));

% 4. Create coordinate grids for the 2D slice (x vs z)
% Transpose P_slice so dimensions match meshgrid(x, z) orientation
[X_2D, Z_2D] = meshgrid(x, z);

% 5. Graph using pcolor
figure;
pcolor(X_2D, Z_2D, P_slice'); 
shading interp; % Smooths the color transitions between grid cells
colorbar;
xlabel('X');
ylabel('Z');
title(['Pressure slice in XZ-plane at Y = ', num2str(y(y_idx))]);
axis equal tight;
