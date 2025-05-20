%Geometry
Ninlet = 60; % points per inlet unit length

xm = -pi;  xM = pi; % range in x
ym = -1;  yM = 1;  % range in y

L = xM-xm;  H = yM-ym; % lengths of x and y domain
hH = H/Ninlet;      % dist between points in y
hL = L/ceil(L/hH); % dist between points in x

% Set up the channel  -L/2<x<L/2, -H/2<y<H/2
x = (xm+hL/2 : hL : xM-hL/2)';  NL = length(x)-1;
y = (ym+hH/2 : hH : yM-hH/2)';  NH = length(y)-1;

[xx, yy] = meshgrid(x,y);

% Darcy number 
Da = 0.3; %Da<1/2 has finitely many real eigvals

[uexact, vexact, pexact] = permeablechannelexact(xx, yy, Da);

% Quiver plot
figure
quiver(xx, yy, uexact, vexact)
axis equal;

% Plotting streamlines with pressure contour
figure
h = contourf(xx, yy, pexact, 12, 'EdgeColor', 'none');
%set(h, 'EdgeColor', 'none');
hold on

stream_in = streamline(xx, yy, uexact, vexact, xx(1:end, 1), yy(1:end, 1));
stream_out = streamline(xx, yy, -uexact, -vexact, xx(1:end, end), yy(1:end, end));

set(stream_in, 'Color', 'black');
set(stream_out, 'Color', 'black');
axis equal;
ax = gca;
ax.FontSize = 14;

colorbar
colormap bone
% FB has logarithmic scale for pressure that gives different contours.
% We think that the discrepancy is that Figure 1 is in dimensional
% variables while we are using the nondimensional. To do: try using
% eq (5) to redimensionalize pressure and see if we get the same
% values as FB.

%%
function [usol, vsol, psol] = permeablechannelexact(x, y, Da)

% solve eq (30) in FB
f=@(lam) tan(lam)^2 - tan(lam)/lam + 1 - 2*Da;
lam1 = fzero(f, sqrt(2*Da));

% set up g_1(y) and derivative, eq (40) in FB
g = cot(lam1)*(Da - 0.5)*sin(lam1*y) + 0.5*y.*cos(lam1*y);
g_prime = cot(lam1)*(Da - 0.5)*lam1*cos(lam1*y) + 0.5*cos(lam1*y) - 0.5*lam1*y.*sin(lam1*y);

% exact solution from eq (49)-(51) in FB
C = -1; % arbitrary constant in solution
psol = C*sinh(lam1*x).*cos(lam1*y);
usol = -(C/lam1)*cosh(lam1*x).*g_prime;
vsol = C*sinh(lam1*x).*g;

end