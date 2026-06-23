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