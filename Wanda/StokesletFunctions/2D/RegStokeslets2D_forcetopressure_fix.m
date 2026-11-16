function [p] = RegStokeslets2D_forcetopressure_fix(y,f,x,ep,blob_num, kappa, normal,ds)

% Computes pressures at set of points when given a set of points and
% forces at those points using the Method of Regularized Stokeslets 
% Based on Cortez, Fluids, 2021
% we are assuming the source points are an ellipse

% Force, normals, b(s) to pressure 
% Wanda Strychalski July 31, 2024

%y = (y1,y2) source points of force
%f = (f1,f2) forces at those source points (not force density)
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
f1 = f(:,1);
f2 = f(:,2); 
x1 = x(:,1); 
x2 = x(:,2); 
n1 = normal(:,1);
n2 = normal(:,2);
fn = dot(f,normal,2); % check the size

inside = isPointInPolygon(x1, x2, y1, y2);

% target points
jumpidx = zeros(size(x));
% [ig,jg] = size(y);
% xf=y(:,1);   yf=y(:,2); %NL=fix(length(yf)*4/11);
% for k1 = 1 : ig
%     for k2 = 1: jg
%         tmp = sqrt((xg(k1,k2)-xf(1:2*NL)).^2+(yg(k1,k2)-yf(1:2*NL)).^2);
%          jk = find( tmp==min(tmp) , 1 );
%          jumpidx(k1,k2) = jk;
%     end
% end


for k1 = 1:M
  [~, ~, minIdx] = findClosestVertex(x1(k1), x2(k1), y1, y2);
  jumpidx(k1) = minIdx;
end
    

%initializing the velocity 
p = zeros(M,1);

% Look over the target points
for j = 1:M
    kappa0 = kappa(jumpidx(j));
    %loop over source points
    for k = 1:N

        %distance between target and source points
        XY1 = x1(j) - y1(k);
        XY2 = x2(j) - y2(k);
        R2 = XY1.^2 + XY2.^2 + ep^2;
        R = sqrt( R2 );

        %computing the pressure
        %Note: S1 = G'(r)/r
        %p(r) = (f . x)S1
        [~, ~, S1, ~, ~] = reg_fncs_withdoublet(ep,R,blob_num);

        fdotXY = ((kappa(k)-kappa0)*n1(k)*XY1 + (kappa(k)-kappa0)*n2(k)*XY2)*ds(k) ;

        p(j) = p(j) +fdotXY.*S1;

    end
    p(j) = p(j)+kappa0*inside(j);
end

   %pg(k1,k2) = sum( fdotx.*Hs/mu )*hL ...
    %           + flag*fn(jumpidx(k1,k2))/1.*(yg(k1,k2)>0).*(yg(k1,k2)<H);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function isPointInPolygon
% px,py: 1D vectors or points (query points)
% xv,yv: 1D vectors containing the coordinates of the 
% curve's vertices.
% Oupute: 1 x length(px) 0 *outside, 1 inside
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function retval = isPointInPolygon(px, py, xv, yv)
% ISPOINTINPOLYGON Determines if a 2D point (px, py) is inside a closed curve.
% xv, yv are vectors containing the coordinates of the curve's vertices.

retval = zeros(length(px),1);

for k=1:length(px(:))

    xt = px(k);
    yt = py(k);

    inside = false;
    numVertices = length(xv);

    % Loop through every edge of the discretized curve
    j = numVertices; % Start with the edge connecting the last vertex to the first
    for i = 1:numVertices
        % Extract vertex coordinates for the current edge
        x1 = xv(i); y1 = yv(i);
        x2 = xv(j); y2 = yv(j);

        % Check if the point's y-coordinate falls within the edge's y-range
        if (y1 > yt) ~= (y2 > yt)
            % Compute the x-coordinate where the horizontal ray intersects the edge
            x_intersect = x1 + (yt - y1) * (x2 - x1) / (y2 - y1);

            % If the intersection point is to the right of the test point, flip state
            if xt < x_intersect
                inside = ~inside;
            end
        end
        j = i; % Move to the next edge
    end
    retval(k) = inside;
end
end
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