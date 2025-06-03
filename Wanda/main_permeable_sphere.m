%main_permeable_sphere

%Based on Example 3 Section 3.1 of Cortez, Fluids, 2021
% Permeable circle with forces = curvature times normal vector
% Compare to an exact solution
% Start with a sphere of radius 1
% Compute the velocity from forces

% Developed for ICERM
% October 1st, 2024
% Author: Wanda Strychalski

close all
clear all
clc


%% Parameters to set

%setting the viscosity
mu = 1;

% initial radius of the sphere
r  = 1;

%%%%%%%%%%%%%%%%%%%
% sphere information
%load md14.0225
load md019.00400
%load md029.00900
%load md28.0841

movie_flag = 0;

if(movie_flag)
    vidObj = VideoWriter('sphere_movie900pts.mp4','MPEG-4');
    open(vidObj);
end

% Update this
pts = md019(:,1:3);
wts = md019(:,4);

X = pts(:,1);
Y = pts(:,2);
Z = pts(:,3);

Tfull = delaunayTriangulation(X, Y, Z);
[T, Xt] = Tfull.freeBoundary();
triN = length(T);
ptN = length(Xt);
N = ptN;

% For computing forces on the sphere
triNormalVec = zeros(triN,3);
triArea = zeros(triN,1);
pointArea = zeros(ptN,1);
[triArea, pointArea] = computeTriangleAreas(Xt,T);

[forces, normalVecs] = computeSphereForces(Xt,T,triArea);

%quiver3(Xt(:,1),Xt(:,2),Xt(:,3),forces(:,1),forces(:,2),forces(:,3));
%axis equal


dA = mean(triArea);
dA = mean(wts);  % Update here!
dA0 = dA;
pointArea = wts;

% write function to compute the curvature and normal vectors
% Blob size: square root of the mean area for the blob size.
%regularization parameter
ep =sqrt(dA);
mu = 1;
%%%%%%%%%%%%%%%%%%%%


blob_num = 1; % This is the only blob where we know the exact solution at this time!

alpha = 1e-1; % was 1e-2

%tfinal = r^2*mu/(2*alpha); % 5, the last time that the solution is valid
X = Xt(:,1);
Y = Xt(:,2);
Z = Xt(:,3);

R = 1/ptN*sum(sqrt(X.^2 + Y.^2+Z.^2));

% source points y
% target points x
y = [X,Y,Z];
y1 = X;
y2 = Y;
y3 = Z;

beta = 4/3*ep*alpha*ones(N,1); % change here! (0.83)
% beta = 1/2*ep*alpha*ones(N,1); % change here!
% b(s)= kappa times beta
% b(s) = beta (f dot n)
%kappa = sqrt(forces(:,1).^2+forces(:,2).^2+forces(:,3).^2);
%[u_beta] = RegStokeslets3D_gforceto_velocity_permeable(y,...
%    y,forces,ep,mu,blob_num, beta, normalVecs);

%u_beta = -u_beta.*pointArea;
%mag_u_beta = sqrt(u_beta(:,1).^2+u_beta(:,2).^2+u_beta(:,3).^2);
%histogram(mag_u_beta);

%vdotn = dot(u_beta,normalVecs); % for checking with the paper
%histogram(vdotn)

R0 = 1;
dt = 1e-3;
tfinal_theory = R0^2*mu/(4*alpha); % 2.5, the last time that the solution is valid

tfinal = tfinal_theory;
out_every = 50;
nt = round(tfinal/dt);
saveNo = floor(nt/out_every)+1;
counter = 1;

R = 1/ptN*sum(sqrt(y1.^2 + y2.^2+y3.^2));

timeArray = 0:dt*out_every:tfinal;
Rexact = sqrt(R0^2-4*alpha/mu*timeArray); % exact solution

% try the first time step with no update using the spherical coordinates
% Convert to spherical coordinate to exactly compute the normal vectors and
% curvature
rho = sqrt(X.^2 + Y.^2 + Z.^2);
theta = atan2(Y,X);
phi = acos(Z./rho);

% check point below
% scatter3(rho.*sin(phi).*cos(theta),rho.*sin(phi).*sin(theta),rho.*cos(phi));

exact_normals = [cos(theta).*sin(phi),sin(theta).*sin(phi),cos(phi)];
% figure;
% quiver3(X,Y,Z,exact_normals(:,1),exact_normals(:,2),exact_normals(:,3))


Rtime = zeros(saveNo,1);
%timeArray = zeros(saveNo,1);
Rtime(counter) = R;
counter = counter+1;

X0 = X;
Y0 = Y;
Z0 = Z;


%looping over time
for j=1:nt % was nt

    time = j*dt;

    % Compute the area of each triangle
    [triArea, pointArea] = computeTriangleAreas(Xt,T);
    dA = R^2*dA0;     % Update here!
    pointArea = R^2*wts;  % Update here!
    %dA = mean(triArea);

    % Compute forces numerically
    % [forces, normalVecs] = computeSphereForces(Xt,T,triArea);
    %forces = (2/R)*normalVecs; % override the forces
    %forces2 = forces.*pointArea;

    % figure(2)
    % quiver3(X,Y,Z,forces(:,1),forces(:,2),forces(:,3))
    % axis equal

    forces = (2/R)*exact_normals; % checking with the exact normals.

    % source points y
    % target points x
    y = [X,Y,Z];
    y1 = X;
    y2 = Y;
    y3 = Z;

    % Compute the porous slip velocity
    [u_beta] = RegStokeslets3D_gforceto_velocity_permeable(y,...
        y,forces,ep,mu,blob_num, beta, normalVecs);

    u_beta = u_beta.*pointArea;
    mag_u_beta = sqrt(u_beta(:,1).^2+u_beta(:,2).^2+u_beta(:,3).^2);
    % should be (alpha/mu)*2*kappa (Laplace's Law + p.6 in the Fluids paper)
    exact_u = (alpha/mu)*2/R;
    comp_u = mean(mag_u_beta);

    % figure(2)
    % quiver3(X,Y,Z,u_beta(:,1),u_beta(:,2),u_beta(:,3))
    % axis equal

    % Cheating section!!
    % Use the normal component of u_beta to update the points
    ubdotn = dot(u_beta,normalVecs,2); % normal component of the velocity
    u_mean = mean(ubdotn); % average of the normal component
    u_beta = u_mean.*normalVecs; %(u_beta.normals)*normals 
    %u_beta = ubdotn.*normalVecs;

    % b = beta.*forces; 
    % b(s) = beta (f dot n) = beta 2 kappa
    %mag_normals = sqrt(normalVecs(:,1).^2+normalVecs(:,2).^2+normalVecs(:,3).^2);

    
    % for checking with the paper
    vdotn = dot(u_beta,normalVecs,2); % for checking with the paper
    kappa = sqrt(forces(:,1).^2+forces(:,2).^2+forces(:,3).^2);
    mag_f = sqrt(forces(:,1).^2+forces(:,2).^2+forces(:,3).^2);

    % Do NOT update when using the exact normals!!
    y1 = y1 + dt*(u_beta(:,1));
    y2 = y2 + dt*(u_beta(:,2));
    y3 = y3 + dt*(u_beta(:,3));

    X = y1;
    Y = y2;
    Z = y3;

    % dynamically update epsilon for accuracy
    R = 1/ptN*sum(sqrt(X.^2 + Y.^2+Z.^2));
    
    ep = 1.5*sqrt(dA);
    beta = 4/3*ep*alpha*ones(N,1);

    if(isnan(R) | R< 1e-6 |R> 1.1)
        fprintf('Stop early\n');
        break;
    end


    if( mod(j,out_every)==0 )
    %fprintf('Computed |u| = %g, Exact |u| = %g\n',exact_u,comp_u);

        Rtime(counter) = R;
        %timeArray(counter) = time;
        counter = counter+1;
        figure(1)
        scatter3(X0,Y0,Z0)
        hold on
        % scatter3(X,Y,Z)
        axis equal
        scatter3(y1,y2,y3,'filled')
        title(sprintf('Time = %4.2f s',time),'fontweight','normal');
        set(gca,'fontsize',16);
        set(gcf,'color','w');
        set(gca,'fontname','Times New Roman'); box on;
        hold off

        if(movie_flag)
            currFrame = getframe(gcf);
            writeVideo(vidObj,currFrame);
        end
        %figure(1)
          %  histogram(mag_u_beta);

         %figure(2)
         % histogram(vdotn)

        percentComplete = time/tfinal*100;
        fprintf('%18.2f%% complete\n',percentComplete);
    end


    Xt = [X Y Z];
end

if(movie_flag)
    close(vidObj);
end

figure;
plot(timeArray,Rexact,'k-');
hold on
grid on
plot(timeArray,Rtime,'bo');
set(gca,'fontsize',16);
set(gcf,'color','w');
set(gca,'fontname','Times New Roman'); box on;
hold off;
legend('Exact solution','Computed');
title([num2str(ptN) ' Points']);
xlabel('Time');
ylabel('Radius');

