function [triArea, pointArea] = computeTriangleAreas(points,triangles)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
Xt = points;
T = triangles;
triN = length(T);
triArea = zeros(triN,1);
ptN = length(Xt);
pointArea = zeros(ptN,1);

totalArea = 0;
for j=1:triN
    ind = T(j,:);
    x0 = Xt(ind(1),:);
    x1 = Xt(ind(2),:);
    x2 = Xt(ind(3),:);

    u = x1-x0;
    v = x2-x1;
    w = x2-x0;

    a = norm(u);
    b = norm(v);
    c = norm(w);
    s = 0.5*(a+b+c);
    area = sqrt(s*(s-a)*(s-b)*(s-c));

    triArea(j) = area;
    pointArea(ind) = pointArea(ind)+1/3*area;
    totalArea = totalArea+area;
end
%diff = 4*pi - totalArea;
%fprintf("Area diff: %f\n", diff);
%tri_avg = mean(triArea);
%fprintf("Triangle area: %6e\n", tri_avg);

end