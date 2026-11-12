%% Plot velocity throughout permeable region
L = 3;
c = 2;
H = 1;
eta = 0.74;
G = 4;
xm = (L+2*c)/2;
Da = 0.3;

x = linspace(0, L+2*c, 100);
y = linspace(0, 2*H, 100);

[xx, yy] = meshgrid(x,y);

c1 = -G/(eta*cosh(eta*L/2)*cos(eta*H));
v_perm = @(x1,x2) c1*sinh(eta*(x1-xm))*((Da - H/2)*cot(eta*H)*sin(eta*(x2-H)) + (x2-H)*cos(eta*(x2-H))/2);

close all
figure
surf(xx,yy,v_perm(xx,yy), 'EdgeColor','none')
view(2)
axis equal

%% Does p(x,y) match p_T at y=0 and y=2H?
params.G = 10;
params.L = 3;
params.c = 2;
params.H = 1;
params.Da = 0.3;
L = params.L;
c = params.c;
G = params.G;
Da = params.Da;

% Newton's method to find lambda1
L1 = sqrt(2*Da);
f = @(L) tan(L)^2-tan(L)/L+1-2*Da;
params.eta = fzero(f,L1);
eta = params.eta;

x1 = linspace(0, L+2*c, 500);
x2 = zeros(size(x1));

p_Fourier = exactPressureFun_Fourier(x1, x2, params);

% Set up for plot
xm = (L+2*c)/2; % midpoint of domain
xL = x1(x1<c & x1>=0);
xC = x1(x1>=c & x1<=L+c);
xR = x1(x1>L+c & x1<=L+2*c);
pL = G*(c-xL+tanh(eta*L/2)/eta);
pC = -G*sinh(eta*(xC-xm))/(eta*cosh(eta*L/2));
pR = G*(L+c-xR-tanh(eta*L/2)/eta);

p_T = [pL, pC, pR];

close all
figure
plot(x1, p_Fourier, '-b', x1, p_T, '--r', 'LineWidth', 2)