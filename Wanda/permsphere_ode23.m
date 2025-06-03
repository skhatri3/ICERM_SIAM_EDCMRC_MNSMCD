%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Developed for ICERM
% Permeable sphere in 3D, where F = kappa*normals
% Code to simulate example 3.1 Fluids (2021) in 3D
% Time update: ODE23
% Date: May 27, 2025
% Author: Wanda Strychalski
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%close all
%clear all
clc

%% Parameters to set

mu = 1;       % fluid viscosity
alpha = 1e-1; % permeability parameter
r  = 1;       % initial radius of the sphere

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sphere information
% Point sets downloaded from: 
% https://web.maths.unsw.edu.au/~rsw/Sphere/Images/MD/md_data.html
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% md14.0225
% md029.00900
% md019.00400
data = load('md029.00900');

% Update this
pts = data(:,1:3); % [x,y,z] points
wts = data(:,4);   % Quadrature weights at each point 
pointArea = wts;   


X = pts(:,1);
Y = pts(:,2);
Z = pts(:,3);
ptN = length(X);
R = 1/ptN*sum(sqrt(X.^2 + Y.^2+Z.^2)); % initial radius

dA = mean(pointArea); % average quadrature weight dA
dA0 = dA;             % initial average quadrature weight
dA0 = 4*pi/ptN;       % try using a uniform weight
pointArea = dA0*ones(ptN,1);% try using a uniform weight

% write function to compute the curvature and normal vectors
% Blob size: square root of the mean area for the blob size.
% regularization parameter
ep = sqrt(dA); % epsilon or delta for MRS

blob_num = 1; % This is the only blob where we know the exact solution at this time!

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Exact solution
%%%%%%%%%%%%%%%%%%%%%%%%%%%
R0 = 1;
dt = 1e-3;
tfinal_theory = R0^2*mu/(4*alpha); % 2.5, the last time that the solution is valid

tfinal = tfinal_theory;
out_every = 50;
nt = round(tfinal/dt);
saveNo = floor(nt/out_every)+1;
counter = 1;

timeArray = 0:dt*out_every:tfinal;
Rexact = sqrt(R0^2-4*alpha/mu*timeArray); % exact solution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Convert to spherical coordinate to exactly 
% compute the normal vectors and curvature
rho = sqrt(X.^2 + Y.^2 + Z.^2);
theta = atan2(Y,X);
phi = acos(Z./rho);

% check point below
% scatter3(rho.*sin(phi).*cos(theta),rho.*sin(phi).*sin(theta),rho.*cos(phi));

% check if pointing outward or inward
exact_normals = [cos(theta).*sin(phi),sin(theta).*sin(phi),cos(phi)];

X0 = X;
Y0 = Y;
Z0 = Z;
z0 = [X0(:);Y0(:);Z0(:)]; %row vector
XX = [X,Y,Z];  % row vector
z = XX(:);  % Stack the dependent variables

beta = 4/3*ep*alpha*ones(ptN,1); % beta depends on alpha

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ODE solver
dt = 0.1;
Tfin = 2.5;
% Set options for the ODE solver
options = odeset('RelTol',1e-10,'AbsTol',1e-9);
tic
normals = exact_normals(:);
pointArea0 = pointArea;   % initial point quadrature weight

smooth_flag = 0;

[t,z]=ode23(@(t,z) velo(t,z,beta,mu,blob_num,dA0,normals,pointArea0,alpha,smooth_flag,Tfin),...
    [(0:dt:Tfin)],z0,options);
toc
%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%
% Graphing the output
Rtime = zeros(length(t),1); % Radius over time

for j = 1:1:length(t)
    N = round(length(z)/3);
    time = t(j);

    X = z(j,1:N);
    Y = z(j,N+1:2*N);
    Z = z(j,2*N+1:3*N);

    figure(1)
    scatter3(X0,Y0,Z0) % initial point positions
    hold on
    axis equal
    scatter3(X,Y,Z,'filled')
    title(sprintf('Time = %4.2f s',time),'fontweight','normal');
    set(gca,'fontsize',16);
    set(gcf,'color','w');
    set(gca,'fontname','Times New Roman'); box on;
    hold off

    Rtime(j) = 1/N*sum(sqrt(X.^2 + Y.^2+Z.^2));

    % if mod(j,8)==1, 
    %     plot([x;x(1)],[y;y(1)],'-','LineWidth',2),axis equal,
    % axis([-1.2 1.2 -0.7 0.7]*1. )
    % title(['t = ',num2str(t(j))]),pause(0.1), hold on
    % xlabel('x'),ylabel('y')
    % end
    %approxim
end

% Plot the radius over time
figure(2);
plot(timeArray,Rexact,'k-','linewidth',2);
hold on
grid on
plot(t,Rtime,'ok','markersize',10,'markersize',10,'linewidth',1,'markerfacecolor','#A7A5A5');
set(gca,'fontsize',16);
set(gcf,'color','w');
set(gca,'fontname','Times New Roman'); box on;
hold off;
legend('Exact solution','Computed');
title([num2str(ptN) ' Points']);
xlabel('Time');
ylabel('Radius');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute the velocity using the source doublets
% Input
%   t = time
%   z = column stacked points (x,y,z)
%   dA0 = initial average quadrature weight for the 
%   unit sphere
%   normals = exact normal vectors using spherical coordinates
%   pointArea0 =  initial quadrature weights for the points
%   on the sphere.
%   alpha = permeability parameter
%   smooth_flag = 1, take the normal component of the velocity and average
%   tfinal = end time
% Output
%   vel = column stacked velocity at points (x,y,z)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function vel=velo(t,z,beta,mu,blob_num,dA0,normals, pointArea0, alpha, smooth_flag, tfinal)

N = round(length(z)/3);
X = z(1:N); 
Y = z(N+1:2*N);
Z = z(2*N+1:3*N);

vel = zeros(3*N,1); % initialize the output 
R = 1/N*sum(sqrt(X.^2 + Y.^2+Z.^2));
kappa_vec = 1./(sqrt(X.^2 + Y.^2+Z.^2)); % curvature at each point as 1/R.

% Updating dA (quadrature weight) as the radius changes
dA = R^2*dA0;   % this is an average scalar
pointArea = R^2*pointArea0; % pointwise quadrature weights

% dynamically update epsilon for MRS for accuracy
ep = 1.5*sqrt(dA); 

if blob_num==1 
    beta = 4/3*ep*alpha*ones(N,1);
end

% Compute the force as 2/R*normals
exact_normals = [normals(1:N), normals(N+1:2*N), normals(2*N+1:3*N)];
%forces = (2/R)*exact_normals; % using average curvature
forces = (2*kappa_vec).*exact_normals; % checking with the exact normals.
%fdotn = dot(forces,normalVecs,2); % normal component of the velocity

y = [X,Y,Z];
% Compute the porous slip velocity
[u_beta] = RegStokeslets3D_gforceto_velocity_permeable(y,...
    y,forces,ep,mu,blob_num, beta, exact_normals);
u_beta = u_beta.*pointArea; % Why do I need an area weighting?
%u_beta_norm = sqrt(u_beta(:,1).^2+u_beta(:,2).^2+u_beta(:,3).^2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Cheating section!!
% The points don't have the same norm unless this section 
% is added in
if(smooth_flag)
    normalVecs = exact_normals;
    % Use the normal component of u_beta to update the points
    ubdotn = dot(u_beta,normalVecs,2); % normal component of the velocity
    %u_beta = ubdotn.*normalVecs;
    u_mean = mean(ubdotn); % average of the normal component
    u_beta = u_mean.*normalVecs; %(u_beta.normals)*normals
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
exact_u = -(alpha/mu)*2/R;
vdotn = dot(u_beta,exact_normals,2); % for checking with the paper
u_error = abs(exact_u-mean(vdotn));
%fprintf('Error: %f\n',u_error);


% plot(abs(ubdotn)-exact_u*ones(N,1))
% %hold on 
% %plot(exact_u*ones(N,1));
% title(sprintf('Time = %4.2f s',t),'fontweight','normal');
% xlabel('Point number');
% set(gca,'fontsize',16);
% set(gcf,'color','w');
% set(gca,'fontname','Times New Roman'); box on;
% legend('$\vec{U} \cdot \vec{n}$','True membrane $\vec{U}$','interpreter','latex')
% hold off
% pause(0.1);


%fprintf('time = %g \n',t);
percentComplete = t/tfinal*100;
%fprintf('%18.2f%% complete\n',percentComplete);

% velocity needs to be Nx1
vel = u_beta(:);

end
