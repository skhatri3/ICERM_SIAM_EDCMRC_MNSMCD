% Plot pressure at the boundaries y=0 and y=2H.

% Parameters
Da = 0.4;
H = 1;
G = 4;
L = 3;
c = 2;

% Use Newton's method to compute lambda
lambda = sqrt(2*Da);
fun = @(L) (Da-H/2)*L*cot(L*H)*cos(L*H) + cos(L*H)/2 - H*L*sin(L*H)/2;
lambda = fzero(fun, lambda);

% Set up for plot
xm = (L+2*c)/2; % midpoint of domain
xL = linspace(0, c, 200);
xC = linspace(c, L+c, 200);
xR = linspace(L+c, L+2*c, 200);
pL = G*(c-xL+tanh(lambda*L/2)/lambda);
pC = -G*sinh(lambda*(xC-xm))/(lambda*cosh(lambda*L/2));
pR = G*(L+c-xR-tanh(lambda*L/2)/lambda);

% Create plot
figure
plot([0,L+2*c],[0,0],'k') % x-axis
hold on
plot([c,c], [max(abs(pL)), -max(abs(pR))], '--k') % left boundary of permeable region
plot([L+c,L+c], [max(abs(pL)), -max(abs(pR))], '--k') % right boundary of permeable region
plot(xL, pL, 'b', xC, pC, 'b', xR, pR, 'b', 'LineWidth',2); % pressure

% plot options
grid on
ltx = 'latex'; % interpreter for labeling plot
ftsz = 14; % font size
xlabel('$x$', 'Interpreter',ltx, 'FontSize', ftsz)
ylabel('$p(x,0)=p(x,2H)$', 'Interpreter',ltx, 'FontSize', ftsz)
text(c,-1,'$c$', 'Interpreter',ltx ,'FontSize', ftsz)
text(L+c,-1,'$L+c$', 'Interpreter',ltx,'FontSize', ftsz)
text(L+2*c,-1,'$L+2c$', 'Interpreter',ltx,'FontSize', ftsz)
text((L+2*c)/2,-1.5,'$\frac{L+2c}{2}$', 'Interpreter',ltx,'FontSize', ftsz-2)
ax = gca; ax.FontSize = ftsz;