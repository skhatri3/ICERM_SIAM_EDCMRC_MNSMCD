%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute the solution of the single differential equation in 3D
% Matches the bottom of p. 6, and Fig. 2 Fluids 2021.
% Date: 01/23/25
% Author: W. Strychalski
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
clc

r0 = 1; % Initial radius
mu = 1; % Viscosity
delta = 0.1; % Some small value
alpha = 0.1; % Permeabiliy
beta = 4/3*delta* alpha; % This comes from the analytic analysis

% Rprime computed from integrating over the sphere
f = @(t,y) -2*beta/(mu*y)*(3/(4*delta)+ (delta^7 + 7*delta^5*y^2-delta^2*(delta^2 + 4*y^2)^(5/2))...
/(4*delta*y^2*(delta^2+4*y^2)^(5/2)));

% Exact solution (computed analytically)
rsol = @(t) sqrt(r0^2-4*alpha/mu*t);
tend = r0^2*mu/(4*alpha); % last possible time
dt = tend/150;
t = [0:dt:tend];

% Computed solution
[t1,rcomp] = ode23(f,[0 tend],r0);

% Graph the exact and computed solutions
plot(t,rsol(t),'k-','LineWidth',1.5);
hold on
plot(t1,rcomp,'ok','markersize',10,'markersize',10,'linewidth',1,'markerfacecolor','#A7A5A5')
legend('Exact','Computed','Location','southwest');
grid on
set(gcf,'color','w');
set(gca,'FontSize',20);
set(gca,'fontname','Times New Roman'); box on;
%saveas(gcf, 'rprime3D_beta', 'epsc');
