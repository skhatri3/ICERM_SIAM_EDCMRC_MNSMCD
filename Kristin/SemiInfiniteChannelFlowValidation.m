% Plots exact solution for semi-infinite channel flow with permeable boundaries
% from Bernardi, F., Chellam, S., Cogan, N. G., & Moore, M. N. J. (2023).
% Stokes flow solutions in infinite and semi‐infinite porous channels.
% Studies in Applied Mathematics, 151(1), 116-140.

%Geometry
Ninlet = 30; % points per inlet unit length

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

[uexact, vexact, pexact] = permeablechannelexact_semiinf(xx, yy, Da, L/yM);

% Quiver plot
figure
quiver(xx, yy, uexact, vexact, 2.5)
axis equal;
ax = gca;
ax.FontSize = 14;

% Plotting streamlines with pressure contour
figure
h = contourf(xx, yy, pexact, 12, 'EdgeColor', 'none');
%set(h, 'EdgeColor', 'none');
hold on

stream_in = streamline(xx, yy, uexact, vexact, xx(1:end, 1), yy(1:end, 1));

set(stream_in, 'Color', 'black');
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
function [usol, vsol, psol] = permeablechannelexact_semiinf(x, y, Da, L_r)

% solve eq (30) in FB
f=@(lam) tan(lam)^2 - tan(lam)/lam + 1 - 2*Da;
lam1 = fzero(f, sqrt(2*Da));

% coefficient from eq (66) of FB
C1 = (-2*L_r)*lam1/sin(lam1);

% set up g_1(y) and derivative, eq (40) in FB
g = cot(lam1)*(Da - 0.5)*sin(lam1*y) + 0.5*y.*cos(lam1*y);
g_prime = cot(lam1)*(Da - 0.5)*lam1*cos(lam1*y) + 0.5*cos(lam1*y) - 0.5*lam1*y.*sin(lam1*y);

% ground-state approximate solution from eq (64)-(68), (74)-(75) in FB
psol = C1*exp(lam1*(x - L_r)).*cos(lam1*y);
usol = C1*(-exp(lam1*(x - L_r))).*g_prime/lam1;
vsol = C1*exp(lam1*(x - L_r)).*g;

end