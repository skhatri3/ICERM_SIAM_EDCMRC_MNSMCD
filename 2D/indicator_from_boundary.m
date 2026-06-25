function chi = indicator_from_boundary(xg, yg, xb, yb)

if xb(1) ~= xb(end) || yb(1) ~= yb(end)
    xb = [xb; xb(1)];
    yb = [yb; yb(1)];
end

inside = inpolygon(xg, yg, xb, yb);

chi = double(inside);

end
