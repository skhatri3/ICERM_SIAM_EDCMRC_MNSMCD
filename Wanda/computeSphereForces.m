function [forces, normalVecs] = computeSphereForces(points,triangles,triArea)
% forces are curvature times area
% Point curvature times normal vector
% This is a force, not a force density

Xt = points;
T = triangles;
triN = length(T);
ptN = length(Xt);

normalVecs = zeros(ptN,3);
hp_vec_data = zeros(ptN,3);
ap_vec_data = zeros(ptN,3);
forces = zeros(ptN,3);

for j=1:triN
    ind = T(j,:);
    x0 = Xt(ind(1),:);
    x1 = Xt(ind(2),:);
    x2 = Xt(ind(3),:);
    u = x1-x0;
    v = x2-x1;
    n1 = cross(u,v);
    n = n1/norm(n1);
    triNormalVec(j,:) = n;
    normalVecs(ind,:) =  normalVecs(ind,:)+1/3*n;
    ap_vec_data(ind,:) = ap_vec_data(ind,:)+1/3*n*triArea(j);

    % curvature at point x0
    zL = x1-x0;
    zR = x2-x0;
    v = zR-zL;
    
    vecH = 0.5*cross(n,v);
    hp_vec_data(ind(1),:) = hp_vec_data(ind(1),:)+vecH;

    % curvature at point x1
    zL = x2-x1;
    zR = x0-x1;
    v = zR-zL;
    
    vecH = 0.5*cross(n,v);
    hp_vec_data(ind(2),:) = hp_vec_data(ind(2),:)+vecH;

    % curvature at point x2
    zL = x0-x2;
    zR = x1-x2;
    v = zR-zL;
    
    vecH = 0.5*cross(n,v);
    hp_vec_data(ind(3),:) = hp_vec_data(ind(3),:)+vecH;

end
%%%%%%%%%%%%%%%%%%%%%%%%%
% compute the curvature
%%%%%%%%%%%%%%%%%%%%%%%%%
for j=1:ptN
    kappa(j) = norm(hp_vec_data(j,:))/norm(ap_vec_data(ind,:));
    lengthN = norm(normalVecs(j,:));
    normalVecs(j,:) = normalVecs(j,:)/lengthN;
    forces(j,:) = 2*kappa(j)*normalVecs(j,:);

end
kappa_avg = mean(kappa);

end