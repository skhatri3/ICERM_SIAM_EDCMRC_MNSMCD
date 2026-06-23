function [p] = RegStokeslets2D_forcetopressure_fix(y,kappa,x,ep,blob_num, b,...
     normal, ds_s)

% Computes pressures at set of points when given a set of points and
% forces at those points using the Method of Regularized Stokeslets 
% Based on Cortez, Fluids, 2021

%Brittany: 
% Changing to have an indicator function input and target points arbitrary

% Force, normals, b(s) to pressure 
% Wanda Strychalski July 31, 2024

%y = (y1,y2) source points of force
%f = (f1,f2) forces density at those source points 
%x = (x1,x2) target points 
%p = pressure evaluated at those target points 
%mu is the viscosity 
%ep is the width of the regularization 
% b(s) is related to the permeability
% normal constains the normal vectors

N = size(y,1); %number of source points 
M = size(x,1); %number of target points


%unpacking the inputs 
y1 = y(:,1);
y2 = y(:,2); 
x1 = x(:,1); 
x2 = x(:,2); 


chi_inside = -indicator_from_boundary(x1, x2, y(:,1), y(:,2));


kappa_s0=zeros(length(x(:,1)),1);
for i=1:length(x(:,1))
    [cx, cy, minIdx] = findClosestVertex(x1(i), x2(i), y1, y2);
    kappa_s0(i)=kappa(minIdx);
end

p=zeros(length(x1), 1);
for i=1:length(x(:,1))
%loop over source points    
for k = 1:N 

    %distance between target and source points 
    XY1 = x1(i) - y1(k); 
    XY2 = x2(i) - y2(k);  
    R2 = XY1.^2 + XY2.^2 + ep^2; 
    R = sqrt( R2 ); 

    %computing the pressure
    %Note: S1 = G'(r)/r
    %p(r) = (f . x)S1
    [~, ~, S1, ~, ~] = reg_fncs_withdoublet(ep,R,blob_num);

    fdotXY = (kappa(k)-kappa_s0(i))*normal(k,1)*XY1 + ...
        (kappa(k)-kappa_s0(i))*normal(k,2)*XY2; 
    p(i) = p(i) + fdotXY.*S1*ds_s(k);
end

end
   
    p=p+kappa_s0.*chi_inside;

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [cx, cy, minIdx] = findClosestVertex(px, py, xv, yv)
% Find the closest discrete vertex in the arrays (xv, yv) to an arbitrary point (px, py)

% Ensure vectors are column vectors
xv = xv(:);
yv = yv(:);

% Compute the squared Euclidean distance to every vertex
distSq = (px - xv).^2 + (py - yv).^2;

% Find the index of the vertex with the minimum distance
[~, minIdx] = min(distSq);

% Return the actual coordinates of that point
cx = xv(minIdx);
cy = yv(minIdx);
end


function chi = indicator_from_boundary(xg, yg, xb, yb)

if xb(1) ~= xb(end) || yb(1) ~= yb(end)
    xb = [xb; xb(1)];
    yb = [yb; yb(1)];
end

inside = inpolygon(xg, yg, xb, yb);

chi = double(inside);

end

