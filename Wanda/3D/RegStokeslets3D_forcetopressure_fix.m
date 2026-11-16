function [p] = RegStokeslets3D_forcetopressure_fix(y, kappa, x, ep, blob_num, normal, dS_s, surface_faces)
% Computes regularized pressure at target 3D points (x) from surface source points (y)
%
% Inputs:
%   y             : N x 3 matrix of source surface points [y1, y2, y3]
%   kappa         : N x 1 array of tension/kappa values at source points
%   x             : M x 3 matrix of target points [x1, x2, x3]
%   ep            : Regularization parameter delta (epsilon)
%   blob_num      : Blob choice index
%   normal        : N x 3 matrix of unit outward normal vectors
%   dS_s          : N x 1 vector of quadrature surface area elements
%   surface_faces : (Optional) F x 3 face triangulation indices for 3D inside check

N = size(y, 1);
M = size(x, 1);

% Unpack target and source points
x1 = x(:,1); x2 = x(:,2); x3 = x(:,3);
y1 = y(:,1); y2 = y(:,2); y3 = y(:,3);

% 1. Determine 3D indicator function chi_inside (-1 inside, 0 outside)
% Replace with your preferred 3D point-in-mesh test (e.g., inpolyhedron or alphaShape)
 %chi_inside = indicator_from_boundary_3D(x, y, surface_faces);
% 1. Determine 3D indicator function chi_inside (-1 inside, 0 outside)
r_sphere = mean(sqrt(y1.^2 + y2.^2 + y3.^2));
%chi_inside = indicator_from_boundary_3D(x, y, surface_faces, r_sphere);
chi_inside = indicator_from_boundary_3D_BEM(x, y, normal, dS_s);

% 2. Find closest vertex kappa_s0 for target points
kappa_s0 = zeros(M, 1);
for i = 1:M
    [~, ~, ~,minIdx] = findClosestVertex3D(x1(i), x2(i), x3(i), y1, y2, y3);
    kappa_s0(i) = kappa(minIdx);
end

% 3. Evaluate 3D regularized pressure integral
p = zeros(M, 1);
for i = 1:M
    p_sum = 0;
    for k = 1:N
        % 3D Euclidean distance + regularization parameter squared
        XY1 = x1(i) - y1(k);
        XY2 = x2(i) - y2(k);
        XY3 = x3(i) - y3(k);
        R2 = XY1^2 + XY2^2 + XY3^2 + ep^2;
        R = sqrt(R2);

        % 3D Regularized pressure gradient kernel S1(r, ep)
        % For standard 3D Cortez blob: S1 = (2*R2 + 3*ep^2) / (4*pi * R2^(5/2))
       [~, ~, S1, ~, ~] = reg_fncs_3D(ep,R,blob_num);

        % Normal directional projection
        fdotXY = -((kappa(k) - kappa_s0(i)) * normal(k,1) * XY1 + ...
                   (kappa(k) - kappa_s0(i)) * normal(k,2) * XY2 + ...
                   (kappa(k) - kappa_s0(i)) * normal(k,3) * XY3);

        p_sum = p_sum + fdotXY * S1 * dS_s(k);
    end
    p(i) = p_sum;
end

% Add jump term for points inside the discretized domain
p = p + kappa_s0 .* chi_inside;

end

%% --- Helper Functions --- %%

% function S1 = reg_fncs_3D_S1(ep, R, blob_num)
% % Returns 3D regularized pressure kernel S1
% % Standard 3D regularized stokeslet pressure kernel derivation:
% r2 = R^2 - ep^2;
% S1 = (2*r2 + 5*ep^2) / (4 * pi * (R^5)); 
% end

function [cx, cy, cz, minIdx] = findClosestVertex3D(px, py, pz, xv, yv, zv)
% Computes 3D squared Euclidean distance to find closest boundary vertex
distSq = (px - xv).^2 + (py - yv).^2 + (pz - zv).^2;
[~, minIdx] = min(distSq);

cx = xv(minIdx);
cy = yv(minIdx);
cz = zv(minIdx);
end

% function chi = indicator_from_boundary_3D(x, y, faces)
% % Evaluates inside/outside indicator using MATLAB alphaShape
% % x: M x 3 target points
% % y: N x 3 surface vertices
% 
% shp = alphaShape(y(:,1), y(:,2), y(:,3));
% inside = inShape(shp, x(:,1), x(:,2), x(:,3));
% 
% chi = -double(inside); % -1 inside, 0 outside
% end

function chi = indicator_from_boundary_3D_BEM(x, y, normal, dS_s)
% Computes 3D indicator function chi(x) via Double-Layer Boundary Integral
% Inputs:
%   x      : M x 3 array of target grid points
%   y      : N x 3 array of surface quadrature points
%   normal : N x 3 matrix of outward unit normal vectors
%   dS_s   : N x 1 vector of quadrature surface area elements
% Output:
%   chi    : M x 1 array (-1 inside, 0 outside)

M = size(x, 1);
chi = zeros(M, 1);

% Use a small regularization distance to prevent singular spikes near boundary nodes
ep_bem = 1e-3; 

for i = 1:M
    % Vector pointing from surface source point y to target field point x
    r_vec = x(i, :) - y;                       % N x 3
    r_sq  = sum(r_vec.^2, 2) + ep_bem^2;       % N x 1 (regularized distance)
    r_dist = sqrt(r_sq);                       % N x 1
    
    % Dot product: (x - y) . normal(y)
    dot_prod = sum(r_vec .* normal, 2);         % N x 1
    
    % Solid angle integrand: ((x - y) . n) / ||x - y||^3
    integrand = dot_prod ./ (r_sq .* r_dist);
    
    % Discrete double-layer surface integral over boundary
    % Evaluates to approx -1 inside (since outward normal n points away from domain),
    % approx -0.5 on boundary, and 0 outside.
    val = (1 / (4 * pi)) * sum(integrand .* dS_s);
    
    % Threshold: -1 inside (val < -0.5), 0 outside
    chi(i) = -double(val < -0.5);
end

end