% inputs: L, c0 (desired c), dx0 (desired spacing)
N = 160;
L = 3;
c0 = 2;
n0 = N*L/(L + 2*c0); % number of spacings between c and L+c
m0 = (N - 1 - n0)/2; % number of grid points before c and after L+c
m = min(round(m0), floor((N-2)/2));
n  = N - 1 - 2*m;  
dx = L/n;

% after computing m, n, and dx (either option above)
p  = 50 - ceil(log2(N*dx));       % power-of-2 scale so that N*dx*2^p is about 2^50
dx = round(dx*2^p)/2^p;           % dyadic dx, changed by a relative ~N*2^-50

c  = (m + 0.5)*dx;                  % exact
L  = n*dx;                               % exact
x  = ((0:N-1) + 0.5)*dx;          % every x(j) exact