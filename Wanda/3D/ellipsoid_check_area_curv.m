%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check the area weights for the T-design points
% on an ellipsoid
%
% Check the accuracy of curvature on an ellipsoid.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


close all
clear all
clc
addpath('../StokesletFunctions/3D/');



%% Parameters to set
point_num  = 3;
factor = 1;
blob_num = 1;
mu = 1;       % fluid viscosity
%alpha = 1e-1; % permeability parameter

% 2. Map to Ellipsoid Geometry
%r = 1;
a = 3; b = 2;
c = 1;
%a = 1; b = 1; c = 1;

% Surface element integrand
integrand = @(theta, phi) sin(theta) .* sqrt( ...
    a^2 * b^2 * cos(theta).^2 + ...
    c^2 * sin(theta).^2 .* (a^2 * sin(phi).^2 + b^2 * cos(phi).^2) );

% Perform high-accuracy adaptive double integration
opts = {'AbsTol', 1e-14, 'RelTol', 1e-14};
true_area = 8 * integral2(integrand, 0, pi/2, 0, pi/2, opts{:});

fprintf('Surface Area (Numerical): %.14f\n', true_area);

r = (a*b*c)^(1/3); % initial radius of the sphere
req = r;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sphere information
% Point sets downloaded from:
% https://web.maths.unsw.edu.au/~rsw/Sphere/EffSphDes/sf.html
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
point_array = zeros(3,1);
for j=1:4
    point_num = j+1; % start at 2, 400 points
    [xf,yf,zf] = TdesignW(point_num);
    pts = [xf(:),yf(:),zf(:)];
    X = pts(:,1);
    Y = pts(:,2);
    Z = pts(:,3);
    ptN = length(X);
    point_array(j) = ptN;
    % Create a perfect triangulation topology based on the uniform sphere
    tri = boundary(X, Y, Z);

    % This is using the ellipse for reference
[idx, ~] = knnsearch(pts, pts, 'K', 25); % this should depend on number of points.

    X_ellip = a * X;
    Y_ellip = b * Y;
    Z_ellip = c * Z;
    P_ellipsoid = [X_ellip, Y_ellip, Z_ellip];
    X0 = X_ellip;
    Y0 = Y_ellip;
    Z0 = Z_ellip;

    dA = zeros(ptN, 1);

    % Compute area of every individual triangle in the mesh
    v1 = P_ellipsoid(tri(:,1), :);
    v2 = P_ellipsoid(tri(:,2), :);
    v3 = P_ellipsoid(tri(:,3), :);

    % Vector cross product handles 3D triangle area calculation
    cross_prod = cross(v2 - v1, v3 - v1, 2);
    tri_areas = 0.5 * sqrt(sum(cross_prod.^2, 2));

    % Distribute 1/3 of each triangle's area to its three corner vertices
    for t = 1:size(tri, 1)
        nodes = tri(t, :);
        dA(nodes) = dA(nodes) + (tri_areas(t) / 3);
    end

    approx_area(j) = sum(dA);
    error_area(j) = abs(true_area - approx_area(j));
    fprintf('Error in area (Numerical): %.14f\n', error_area);

    z = [X0,Y0,Z0];
    for i = 1:ptN
        target_pt = z(i, :);
        local_neighbors = z(idx(i, :), :);

        [K(i), H(i), Normals(i,:)] = computePointCurvature(target_pt, local_neighbors);
    end
    H = abs(H); % adjust the mean curvature
    H = H(:);

    % --- IMPLICIT ANALYTICAL MEAN CURVATURE ---
    % For an ellipsoid defined by F(x,y,z) = x^2/a^2 + y^2/b^2 + z^2/c^2 - 1 = 0
    % The components of the gradient vector:
    vX = X_ellip / a^2;
    vY = Y_ellip / b^2;
    vZ = Z_ellip / c^2;

    % Length of the gradient vector
    grad_norm = sqrt(vX.^2 + vY.^2 + vZ.^2);

    % Exact Mean Curvature using implicit differential geometry
    num_H = (vX.^2.*(1/b^2 + 1/c^2) + vY.^2.*(1/a^2 + 1/c^2) + vZ.^2.*(1/a^2 + 1/b^2));
    true_H = num_H ./ (2 * grad_norm.^3);

    error2 = abs(H - true_H);
    error_curv(j) = max(error2);
    fprintf('Mean Curvature Error: %e\n', mean(error2));

    % ==========================================
% FIGURE 1: NUMERICAL ESTIMATION (LOCAL QUADRIC)
% ==========================================

figure(1); % Focuses/reuses Figure 1
clf;       % Clear existing plots
set(gcf, 'Name', 'Numerical Mean Curvature', 'Color', 'w');
% Using abs handles the local sign ambiguity of PCA normals
hSurface = trisurf(tri, X_ellip, Y_ellip, Z_ellip, H);

shading interp; lighting gouraud; camlight head; axis equal; grid on;
colormap(jet);        
c1 = colorbar;
c1.Label.String = 'Computed Numerical H';
xlabel('X'); ylabel('Y'); zlabel('Z');
%title('Numerical Local Quadric Fit');
set(gca,'fontname','Times New Roman');
set(gca,'fontsize',16);
%axis([-1.25 1.25 -1.25 1.25 -1.25 1.25]);
%clim([0.151, 0.903]); % Lock color limits to true mathematical bounds
pause(0.1);
end
saveas(gcf,'ellipse_curvature','epsc');
x_values = 1./sqrt(point_array);
x_values = 1./point_array;

figure;
h1 = plot(x_values,error_area,'ko--','LineWidth',1.5,'MarkerFaceColor','k');
hold on
p = polyfit(x_values,error_area,1);
points_refined = linspace(x_values(1),x_values(end),10);
y1 = polyval(p,points_refined);

h2 = plot(points_refined,y1,LineWidth=2);
xlabel('$1/N \sim \ell^2$', 'Interpreter', 'latex');
ylabel('Error');
%title('Error in Ellipsoid Area');
set(gcf,'color','w');
set(gca,'fontname','Times New Roman');
set(gca,'fontsize',16);
grid on
% 4. Format fit equation string (e.g., y = mx + b)
fit_str = sprintf('Linear Fit: $y = %.1f \\cdot 10^{2} x - %.1f \\cdot 10^{-4}$', ...
    p(1)/100, abs(p(2))*1e4);

% 5. Add legend using explicit handles
legend([h1, h2], {'Data Points', fit_str}, ...
       'Interpreter', 'latex', 'Location', 'southeast');
hold off
saveas(gcf,'ellipse_area_converge','epsc');

figure;
h1 = plot(x_values,error_curv,'ko--','LineWidth',1.5,'MarkerFaceColor','k');
hold on
p = polyfit(x_values,error_curv,1);
points_refined = linspace(x_values(1),x_values(end),10);
y1 = polyval(p,points_refined);
h2 = plot(points_refined,y1,LineWidth=2);
xlabel('$1/N \sim \ell^2$', 'Interpreter', 'latex');
ylabel('Error');
%title('Error in Ellipsoid Curvature');
set(gcf,'color','w');
set(gca,'fontname','Times New Roman');
set(gca,'fontsize',16);
grid on
% 4. Format fit equation string (e.g., y = mx + b)
fit_str = sprintf('Linear Fit: $y = %2.0f x - %.1f \\cdot 10^{-2}$', ...
    p(1), abs(p(2))*1e2);

% 5. Add legend using explicit handles
legend([h1, h2], {'Data Points', fit_str}, ...
       'Interpreter', 'latex', 'Location', 'southeast');
hold off
hold off
saveas(gcf,'ellipse_curvature_converge','epsc');


