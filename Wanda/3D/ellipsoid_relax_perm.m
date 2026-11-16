%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Developed for ICERM
% Permeable ellipsoid in 3D, where F = (Req/R(t)-1)*normals
% Code to simulate example 3.2 Fluids (2021) in 3D
% Time update: Forward Euler
% Date: June 12, 2025
% Author: Wanda Strychalski
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% point_num = 1,2,3 or 4 (select point number)
%                         see TdesignW.m
% blob_num = 1 or 3 (1 is less accurate)
% delta_update  = 0, do not update with the radius
%               = 1, scale delta with the current radius
% factor        = ep = 2*factor*sqrt(dA);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% errorData = [t,errorOver time]
% exactData = [timeArray, exact solution]
% computedData = [t, compute solution with
%                source doublets]
% ptN = number of points on the unit sphere
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [errorData,exactData,computedData,ptN] = ...
%     permsphere_eq_ode23(point_num, blob_num, delta_update,factor)
close all
clear all
clc
addpath('../StokesletFunctions/3D/');

video_flag = 0;
if(video_flag)
    vidObj = VideoWriter('ellipsoid_perm_relax_movie.mp4','MPEG-4');
    open(vidObj);
end

%% Parameters to set
point_num  = 3;
factor = 1;
blob_num = 1;
mu = 1;       % fluid viscosity
alpha = 1e-1; % permeability parameter

% 2. Map to Ellipsoid Geometry
r = 1;
a = 1.6; b = 0.8; 
c = r^3/(a*b); 
% a = 1; b = 1; c = 1;  % Sphere

r = (a*b*c)^(1/3); % initial radius of the sphere
req = r;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sphere information
% Point sets downloaded from:
% https://web.maths.unsw.edu.au/~rsw/Sphere/EffSphDes/sf.html
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[xf,yf,zf] = TdesignW(point_num);
pts = [xf(:),yf(:),zf(:)];

% Update this
%pts = data(:,1:3); % [x,y,z] points
%wts = data(:,4);   % Quadrature weights at each point
wts = 0;
pointArea = wts;


X = pts(:,1);
Y = pts(:,2);
Z = pts(:,3);
ptN = length(X);
R = 1/ptN*sum(sqrt(X.^2 + Y.^2+Z.^2)); % initial radius

% Create a perfect triangulation topology based on the uniform sphere
tri = boundary(X, Y, Z);

dA0 = 4*pi/ptN;       % try using a uniform weight
pointArea = dA0*ones(ptN,1);% try using a uniform weight

% write function to compute the curvature and normal vectors
% Blob size: square root of the mean area for the blob size.
% regularization parameter
% I need the 2 here so that the epsilon covers the average distance between
% 2 points on the sphere
ep = 2*factor*sqrt(dA0); % epsilon or delta for MRS

% if(blob_num == 3)
%     ep = ep/2.1545; % to match the height of blob 1
% end

% blob_num = 1; % Choose 1 or 3.

X_ellip = a * X;
Y_ellip = b * Y;
Z_ellip = c * Z;
P_ellipsoid = [X_ellip, Y_ellip, Z_ellip];
X0 = X_ellip;
Y0 = Y_ellip;
Z0 = Z_ellip;
% Assumes P_current is the N x 3 matrix of moved coordinates at the current timestep
% Assumes 'tri' is the static triangulation matrix from the sphere initialization

dA = zeros(ptN, 1);

% Compute area of every individual triangle in the mesh
v1 = P_ellipsoid(tri(:,1), :);
v2 = P_ellipsoid(tri(:,2), :);
v3 = P_ellipsoid(tri(:,3), :);

% Vector cross product handles 3D triangle area calculation
cross_prod = cross(v2 - v1, v3 - v1, 2);
tri_areas = 0.5 * sqrt(sum(cross_prod.^2, 2));

% Distribute 1/3 of each triangle's area to its three corner vertices
for t = 1:size(tri, 1)
    nodes = tri(t, :);
    dA(nodes) = dA(nodes) + (tri_areas(t) / 3);
end

% Verification: sum(w_fsi) will always equal the exact current surface area of the deformed shape


dt = 1e-2;

tfinal = 1.5;
out_every = 10;
nt = round(tfinal/dt);
saveNo = floor(nt/out_every)+1;
counter = 1;

% check point below
% scatter3(rho.*sin(phi).*cos(theta),rho.*sin(phi).*sin(theta),rho.*cos(phi));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ODE solver
%dt = tfinal/100;\
timeArray  = [0:dt:tfinal];

% Set options for the ODE solver
%options = odeset('RelTol',1e-10,'AbsTol',1e-9);
tic

pointArea0 = pointArea;   % initial point quadrature weight

% This is using the ellipse for reference
[idx, ~] = knnsearch(pts, pts, 'K', 25); % this should depend on number of points.

z0 = [X(:);Y(:);Z(:)]; %row vector
z = [X_ellip(:), Y_ellip(:), Z_ellip(:)];

z1 = z(:,1);
z2 = z(:,2);
z3 = z(:,3);

% Updating dA (quadrature weight) as the radius changes
%dA = req^2*dA0;   % this is an average scalar

% dynamically update epsilon for MRS for accuracy
dA0  = mean(dA);
ep = 2*factor*r^2*sqrt(dA0); % Not sure if this is correct

% For the permeable code
if blob_num==1 
    beta = 4/3*ep*alpha*ones(ptN,1);
elseif blob_num == 3
    % ep = ep/2.1545; % to match the height of blob 1
    beta = 4/9*ep*alpha*ones(ptN,1);
else
    beta = 0;
end


% Graphing the output
Rtime = zeros(nt,1); % Radius over time
[max_val, max_idx] = max(z1);
Rtime(1) = z1(max_idx);
% looping over time

[xi_q, eta_q, w_q] = getTriangleQuadrature();
Nq = length(w_q);

plane_pt = [0 0 0];
plane_n = [0 1 0];
segments0 = sliceTrisurf({tri,z1,z2,z3}, plane_pt,plane_n); %for plane y=0.

for j=1:nt

    time = timeArray(j+1);
    %vel=velo(t,z,mu,idx,blob_num,dA0,req,pointArea0,factor, Tfin);
    z = [z1, z2, z3];

    for i = 1:ptN
        target_pt = z(i, :);
        local_neighbors = z(idx(i, :), :);

        [K(i), H(i), Normals(i,:)] = computePointCurvature(target_pt, local_neighbors);
    end
    H = abs(H); % adjust the mean curvature
    H = H(:);
        
    forces = -2*H.*Normals; % checking with the exact normals.
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Initialize your velocity array to zero
    % This is for Gaussian quadrature
    % u = zeros(ptN, 3); 
    % M_tri = size(tri, 1);
    
   
    % 3. Loop over every triangle element
    % for t = 1:M_tri
    %     idx_nodes = tri(t, :);
    % 
    %     % Node positions for this triangle
    %     v1 = z(idx_nodes(1), :);
    %     v2 = z(idx_nodes(2), :);
    %     v3 = z(idx_nodes(3), :);
    % 
    %     % Node force densities for this triangle
    %     f1 = forces(idx_nodes(1), :);
    %     f2 = forces(idx_nodes(2), :);
    %     f3 = forces(idx_nodes(3), :);
    % 
    %     % Area Jacobian: ||(v2-v1) x (v3-v1)||
    %     J = norm(cross(v2 - v1, v3 - v1)); 
    % 
    %     % Loop over the 7 internal Gaussian quadrature points
    %     for q = 1:Nq
    %         % Linear shape functions
    %         N1 = 1 - xi_q(q) - eta_q(q);
    %         N2 = xi_q(q);
    %         N3 = eta_q(q);
    % 
    %         % Interpolate 3D coordinates and force density at the quad point
    %         xq = N1*v1 + N2*v2 + N3*v3;
    %         fq = N1*f1 + N2*f2 + N3*f3;
    % 
    %         % The exact differential area weight for this point force
    %         point_force = fq * (w_q(q) * J);
    % 
    %         % 4. Compute regularized Stokeslet contribution to all target nodes
    %         for i = 1:ptN
    %             target_pt = z(i, :);
    % 
    %             % Distance calculation
    %             XY1 = target_pt(1) - xq(1);
    %             XY2 = target_pt(2) - xq(2);
    %             XY3 = target_pt(3) - xq(3);
    %             R2 = XY1^2 + XY2^2 + XY3^2 + ep^2;
    %             R = sqrt(R2);
    % 
    %             % Evaluate the Cortez regularization functions directly
    %             [A, B] = reg_fncs_3D(ep, R, blob_num);
    %             fdotXY = point_force(1)*XY1 + point_force(2)*XY2 + point_force(3)*XY3;
    % 
    %             % Accumulate the velocity components
    %             u(i, 1) = u(i, 1) + point_force(1)*A + fdotXY*B*XY1;
    %             u(i, 2) = u(i, 2) + point_force(2)*A + fdotXY*B*XY2;
    %             u(i, 3) = u(i, 3) + point_force(3)*A + fdotXY*B*XY3;
    %         end
    %     end
    % end
    % 
    % % Apply fluid viscosity scaling
    % u = u / mu;

    % end Gaussian quadrature
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Compute the force as 2/R*normals
    % H = ones(ptN,1);
    % forces = -(2/R)*exact_normals; % using average curvature

    force_density = forces.*dA;

    force_int = [sum(force_density(:,1)),sum(force_density), sum(force_density(:,3))]; % need a dA
    %fdotn = dot(forces,normalVecs,2); % normal component of the velocity
    %  RegStokeslets3D_forcetovelocity(y,f,x,ep,mu,blob_num)
    u = RegStokeslets3D_forcetovelocity(z,force_density,z,ep,mu,blob_num);
    udotn = dot(u,Normals,2);
    magu = sqrt(u(:,1).^2+u(:,2).^2+u(:,3).^2);
    magf = sqrt(forces(:,1).^2+forces(:,2).^2+forces(:,3).^2);

    [u_beta] = RegStokeslets3D_gforceto_velocity_permeable(z,...
        z,-force_density,ep,mu,blob_num, beta, Normals);

    %u_beta = zeros(size(u));
    u = u+u_beta;

    z1 = z1+dt*u(:,1);
    z2 = z2+dt*u(:,2);
    z3 = z3+dt*u(:,3);
    Rtime(j+1) = z1(max_idx); % largest deformation in x-direction.

    if(mod(j,out_every)==0)

        figure(1)
        scatter3(X0,Y0,Z0) % initial point positions
        hold on
        axis equal
        scatter3(z1,z2,z3,'filled')
        f = forces;
        quiver3(z1,z2,z3,u(:,1),u(:,2),u(:,3),'-k');
        %length = norm()
        hold off
        title(sprintf('Time = %4.4f s',time),'fontweight','normal');
        set(gca,'fontsize',16);
        set(gcf,'color','w');
        set(gca,'fontname','Times New Roman'); box on;
        hold off

        % ==========================================
        % FIGURE 2
        % ==========================================
        figure(2);
        set(gcf, 'Name', 'Numerical u dot n', 'Color', 'w');
        % Using abs handles the local sign ambiguity of PCA normals
        hSurface = trisurf(tri, z1, z2, z3, udotn);
        % Force the bounding box limits to stay identical across all timesteps
        xlim([-6, 6]);
        ylim([-6, 6]);
        zlim([-6, 6]);
        % xlabel('X Axis', 'FontName', 'Times New Roman', 'FontSize', 16);
        % ylabel('Y Axis', 'FontName', 'Times New Roman', 'FontSize', 16);
        % zlabel('Z Axis', 'FontName', 'Times New Roman', 'FontSize', 16);
        % Set the title with a custom FontName and specific FontSize of 16
        title(sprintf('Time = %4.2f s', time), 'FontName', 'Times New Roman', 'FontSize', 16, 'FontWeight', 'normal');
        set(gca,'fontname','Times New Roman');  
        set(gca,'fontsize',16);
        shading interp; lighting gouraud; camlight head; axis equal; grid on;
        colormap(jet);
        c1 = colorbar;
        clim([-0.07,0.07]);
        c1.Label.Interpreter = 'latex';
        c1.Label.String = 'Computed $\mathbf{u} \cdot \mathbf{n}$';
        % Set the interpreter property to latex

        % Optional: Increase the font size for better readability
        c1.Label.FontSize = 14;
        xlabel('X'); ylabel('Y'); zlabel('Z');
        hold off;
        if(video_flag)
            currFrame = getframe(gcf);
            writeVideo(vidObj,currFrame);
        end
        plane_pt = [0 0 0];
        plane_n = [0 1 0];
        figure(3);
        clf
        segments = sliceTrisurf({tri,z1,z2,z3}, plane_pt,plane_n); %for plane y=0.
        figure(4)
        clf
        style = 'ro-';
        plot2DSlice(segments, plane_pt, plane_n,style);
        hold on
        style = '-k';
        plot2DSlice(segments0, plane_pt, plane_n,style);
        axis([-1.2 1.2 -1.2 1.2])
        axis equal
        hold off
        pause(0.1)
    end
  

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    dA = zeros(ptN, 1); % update dA

    z = [z1, z2, z3];
    % Compute area of every individual triangle in the mesh
    v1 = z(tri(:,1), :);
    v2 = z(tri(:,2), :);
    v3 = z(tri(:,3), :);

    % Vector cross product handles 3D triangle area calculation
    cross_prod = cross(v2 - v1, v3 - v1, 2);
    tri_areas = 0.5 * sqrt(sum(cross_prod.^2, 2));

    % Distribute 1/3 of each triangle's area to its three corner vertices
    for t = 1:size(tri, 1)
        nodes = tri(t, :);
        dA(nodes) = dA(nodes) + (tri_areas(t) / 3);
    end

    % dynamically update epsilon for MRS for accuracy
    dA0  = mean(dA);
    R = 1/ptN*sum(sqrt(z1.^2 + z2.^2+z3.^2)); % initial radius
    r = mean(R);
    ep = 2*factor*r^2*sqrt(dA0); % Not sure if this is correct

    % For the permeable code
    if blob_num==1
        beta = 4/3*ep*alpha*ones(ptN,1);
    elseif blob_num == 3
        % ep = ep/2.1545; % to match the height of blob 1
        beta = 4/9*ep*alpha*ones(ptN,1);
    else
        beta = 0;
    end

end

% figure(2);
% plot(timeArray,Rtime,'ok','markersize',10,'markersize',10,'linewidth',1,'markerfacecolor','#A7A5A5');
% grid on
% set(gca,'fontsize',16);
% set(gcf,'color','w');
% set(gca,'fontname','Times New Roman'); box on;
% hold off;
% %legend('Exact solution','Computed');
% title([num2str(ptN) ' Points']);
% xlabel('Time');
% ylabel('Z(t) at point 1');

% [t,z]=ode23(@(t,z) velo(t,z,mu,idx,blob_num,dA0,req,pointArea0,factor, Tfin),...,
%     [(0:dt:Tfin)],z0,options);
 toc

 if(video_flag)
     close(vidObj);
 end

 figure;
 plot(timeArray,Rtime,'k-','LineWidth',1.5);
 grid on
 xlabel('Time');
 ylabel('Max in x-direction');
 set(gca,'fontsize',16);
 set(gcf,'color','w');
 set(gca,'fontname','Times New Roman'); box on;
%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%

% Output data for the function
%computedData = [t,Rtime];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute the velocity using the source doublets
% Input
%   t = time
%   z = column stacked points (x,y,z)
%   idx = neighbors for each point.
%   blob_num = 1,3
%   dA0 = initial average quadrature weight for the 
%   req = resting radius
%   pointArea0 =  initial quadrature weights for the points
%   on the sphere.
%   factor refers to:   ep = 2*factor*sqrt(dA);
%   tfinal = end time
% Output
%   vel = column stacked velocity at points (x,y,z)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function vel=velo(t,z,mu,idx, blob_num,dA0, req, pointArea0, factor, tfinal)

N = round(length(z)/3);
X = z(1:N); 
Y = z(N+1:2*N);
Z = z(2*N+1:3*N);

vel = zeros(3*N,1); % initialize the output 

% Updating dA (quadrature weight) as the radius changes
dA = req^2*dA0;   % this is an average scalar

% dynamically update epsilon for MRS for accuracy
ep = 2*factor*sqrt(dA); % Not sure if this is correct
% if(delta_update == 1)
%     ep = 2*factor*sqrt(dA); % .5*sqrt(dA)
% end

% if blob_num==1 
%     beta = 4/3*ep*alpha*ones(N,1);
% elseif blob_num == 3
%     % ep = ep/2.1545; % to match the height of blob 1
%     beta = 4/9*ep*alpha*ones(N,1);
% else
%     beta = 0;
% end

P_ellipsoid = [X, Y, Z];
for i = 1:size(P_ellipsoid, 1)
    target_pt = P_ellipsoid(i, :);
    local_neighbors = P_ellipsoid(idx(i, :), :);
    
    [K(i), H(i), Normals(i,:)] = computePointCurvature(target_pt, local_neighbors);
end
H = abs(H); % adjust the mean curvature
H = H(:);

% Compute the force as 2/R*normals
%forces = (2/R)*exact_normals; % using average curvature
forces = -2*H.*Normals; % checking with the exact normals.
%fdotn = dot(forces,normalVecs,2); % normal component of the velocity

y = [X,Y,Z];
f = forces;

u = RegStokeslets3D_forcetovelocity(y,f,y,ep,mu,blob_num);

% Compute the porous slip velocity
%[u_beta] = RegStokeslets3D_gforceto_velocity_permeable(y,...
 %   y,forces,ep,mu,blob_num, beta, exact_normals);
%u_beta = u_beta.*pointArea; % Why do I need an area weighting?
%u_beta_norm = sqrt(u_beta(:,1).^2+u_beta(:,2).^2+u_beta(:,3).^2);

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
vel = u(:);

end

function f = rexact(t,alpha,mu,r0,req)
retval = 0;
c = (r0/req-1)*exp(r0/req-1);
retval = req*(1+lambertw(c*exp(-alpha*t/(mu*req))));

f = retval;
end

function segments = sliceTrisurf(Tri, P0, n)
% Tri  - triangulation object or [T, X, Y, Z] where T is Nx3 and X,Y,Z are vectors
% P0   - a point on the slicing plane [x y z]
% n    - plane normal [nx ny nz]
% segments - Mx2x3 array of segment endpoints (M segments)

% Accept different inputs
if isa(Tri,'triangulation')
    TR = Tri;
else
    T = Tri{1}; X = Tri{2}; Y = Tri{3}; Z = Tri{4};
    TR = triangulation(T, [X(:) Y(:) Z(:)]);
end

V = TR.Points;
Faces = TR.ConnectivityList;
P0 = P0(:).';
n = n(:)./norm(n);

segments = [];
for i=1:size(Faces,1)
    idx = Faces(i,:);
    verts = V(idx,:);            % 3x3
    % signed distances to plane
    d = (verts - P0) * n;
    pts = [];
    % check each edge
    edges = [1 2; 2 3; 3 1];
    for e=1:3
        a = edges(e,1); b = edges(e,2);
        da = d(a); db = d(b);
        if da==0
            pts(end+1,:) = verts(a,:); %#ok<AGROW>
        end
        if da*db < 0     % sign change => intersection
            t = da/(da-db);         % linear interp parameter from a to b
            ip = verts(a,:) + t*(verts(b,:)-verts(a,:));
            pts(end+1,:) = ip; %#ok<AGROW>
        elseif db==0
            pts(end+1,:) = verts(b,:); %#ok<AGROW>
        end
    end
    % keep unique intersection points
    if ~isempty(pts)
        pts = unique(round(pts,12),'rows'); % tolerance
        if size(pts,1)==2
            segments(end+1,:,:) = pts; %#ok<AGROW>
        end
    end
end

% Plot trisurf and slice segments if called without outputs
%if nargout==0
    figure(3); hold on; trisurf(TR,'FaceAlpha',0.6,'EdgeColor','k');
    for k=1:size(segments,1)
        seg = squeeze(segments(k,:,:));
        plot3(seg(:,1), seg(:,2), seg(:,3), 'r-', 'LineWidth',2);
    end
    axis equal; view(3); hold off;
%end
end

function plot2DSlice(segments, P0, n,style)
% segments : Mx2x3 array of 3-D segment endpoints
% P0       : 1x3 point on plane
% n        : 1x3 plane normal

P0 = P0(:).';
n = n(:).' / norm(n);

% pick tangent axis u (not parallel to n)
if abs(dot(n, [1 0 0])) < 0.9
    u = cross(n, [1 0 0]);
else
    u = cross(n, [0 1 0]);
end
u = u / norm(u);
v = cross(n, u);        % orthonormal v

% build 2x3 projection matrix (rows are u and v)
R = [u; v];             % 2x3

M = size(segments,1);
XY = zeros(M,2,2);      % Mx2x2: [p1; p2] each as [x y]

for k = 1:M
    pts3 = squeeze(segments(k,:,:));    % 2x3
    pts3 = pts3 - P0;                   % shift origin
    pts2 = (R * pts3.').';              % 2x2
    XY(k,:,:) = pts2;
end

% Plot all segments in 2D
figure(4); hold on; axis equal;
for k = 1:M
    seg = squeeze(XY(k,:,:));          % 2x2
    plot(seg(:,1), seg(:,2), style, 'LineWidth', 1.5);
    %% 
end
xlabel('u'); ylabel('v'); grid on; hold off;   set(gca,'FontSize',20);
set(gcf,'color','w');
set(gca,'fontname','Times New Roman'); box on;
end