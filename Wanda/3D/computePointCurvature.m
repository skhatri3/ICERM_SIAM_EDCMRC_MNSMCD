function [K, H, normal] = computePointCurvature(pt, neighbors)
    % 1. Center neighbors at the origin
    centered = neighbors - pt;
    shape_center = 0;
    
    % 2. Get PCA normal for local coordinate system
    [coeff, ~] = pca(centered);
    z_axis = coeff(:,3);
    x_axis = coeff(:,1);
    y_axis = coeff(:,2);
    
    % 3. Transform neighbors to local (u, v, w) coordinates
    % Local coordinates: u = x-axis, v = y-axis, w = normal-axis
    local_pts = centered * [x_axis, y_axis, z_axis];
    u = local_pts(:,1);
    v = local_pts(:,2);
    w = local_pts(:,3);
    
    % 4. Fit Quadric: w = a*u^2 + b*u*v + c*v^2 + d*u + e*v
    % Construct design matrix
    M = [u.^2, u.*v, v.^2, u, v];
    params = M \ w; % Solve linear least squares
    a = params(1); b = params(2); c = params(3); d = params(4); e = params(5);
    
    % 5. Curvities at the origin (u=0, v=0)
    % The normal vector in local space is (-d, -e, 1) normalized
    local_n = [-d; -e; 1] / sqrt(d^2 + e^2 + 1);
    normal = [x_axis, y_axis, z_axis] * local_n; % Back to global space

    % ----------------------------------------------------
    % FIX: Force normal to point OUTWARD
    % ----------------------------------------------------
    % Vector pointing from the shape's center to the current point
    outward_direction = pt - shape_center; 
    
    % If the dot product is negative, they point in opposite directions
    if dot(normal, outward_direction) < 0
        normal = -normal;
    end
    
    % Fundamental forms for curvature
    E = 1 + d^2; F = d*e; G = 1 + e^2;
    L = 2*a; M_curv = b; N = 2*c;
    
    % Gaussian Curvature (K) and Mean Curvature (H)
    %K = (L*N - M_curv^2) / (E*G - F^2);
    %H = (E*N + G*L - 2*F*M_curv) / (2 * (E*G - F^2));

    % FIX: Incorporate the first-derivative slope terms (d and e)
    denom_K = (1 + d^2 + e^2)^2;
    denom_H = (1 + d^2 + e^2)^1.5;

    % Corrected Curvature Estimations
    K = (4 * a * c - b^2) / denom_K;
    H = (a * (1 + e^2) - b * d * e + c * (1 + d^2)) / denom_H;
end