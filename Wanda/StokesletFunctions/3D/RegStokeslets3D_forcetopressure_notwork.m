function [p] = RegStokeslets3D_forcetopressure(y,f,x,ep,blob_num)

% Computes pressures at set of points when given a set of points and
% forces at those points using the Method of Regularized Stokeslets 
% Based on Cortez, SIAM J. Sci Comput. 2001

% Force to pressure code developed by Wanda Strychalski
% September 20206

%y = (y1,y2) source points of force
%f = (f1,f2) forces at those source points (not force density)
%x = (x1,x2) target points 
%p = pressure evaluated at those target points 
%mu is the viscosity 
%ep is the width of the regularization 

N = size(y,1); %number of source points 
M = size(x,1); %number of target points 

%unpacking the inputs 
y1 = y(:,1);
y2 = y(:,2);
y3 = y(:,3);

f1 = f(:,1);
f2 = f(:,2); 
f3 = f(:,3); 

x1 = x(:,1); 
x2 = x(:,2); 
x3 = x(:,3); 

%initializing the pressure 
p = zeros(M,1);

%loop over source points    
for k = 1:N 

    %distance between target and source points 
    XY1 = x1(:) - y1(k); 
    XY2 = x2(:) - y2(k); 
    XY3 = x3(:) - y3(k);
    
    R2 = XY1.^2 + XY2.^2 + XY3.^2 + ep^2; 
    R = sqrt( R2 ); 

    %computing the pressure
    %Note: S1 = G'(r)/r
    %p(r) = (f . x)S1
    [~, ~, S1, ~, ~] = reg_fncs3D(ep,R,blob_num);

    fdotXY = f1(k)*XY1 + f2(k)*XY2 + f3(k)*XY3; 
    
    p(:) = p(:) + fdotXY.*S1;

end
end
function [cx, cy, minIdx] = findClosestVertex(px, py, pz, xv, yv, zv)
% Find the closest discrete vertex in the arrays (xv, yv) to an arbitrary point (px, py)

% Ensure vectors are column vectors
xv = xv(:);
yv = yv(:);
zv = zv(:);

% Compute the squared Euclidean distance to every vertex
distSq = (px - xv).^2 + (py - yv).^2+(pz - zv).^2;

% Find the index of the vertex with the minimum distance
[~, minIdx] = min(distSq);

% Return the actual coordinates of that point
cx = xv(minIdx);
cy = yv(minIdx);
cz = zv(minIdx);
end