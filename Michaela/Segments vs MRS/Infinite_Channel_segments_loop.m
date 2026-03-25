% Flow through permeable channel using exact solutions on left and right
% boundaries from Bernardi et al 2023.

% This code implements the Segment method using nested for-loops and was
% used to validate the segment force-to-velocity formulations.
%
% Developed by Michaela Kubacki and Ricardo Cortez March 2026


clear all

%% Parameters to set

% Model Parameters
Da = 0.4; % Darcy Number
mu = 1; % Viscosity
blob_num = 2; % blob choice

% Number of segment and target points
%N = 160; % Number of source points (along top  and bottom boundaries)

Ny1 = 148; % Number of y-points on top/bottom boundary (no corners)
% Ny1 + 1 segments on top/bottom
Ny2 = 100;  % Number of y-points on inlet/outlet (includes corners)
% Ny2 - 1 segments on left/right
Nseg = 2*(Ny1+Ny2); % Total number of segments on the boundary

Nx1 = Ny1*0.5; % Number of target points in x direction for full channel calc
Nx2 = Ny2*0.5; % Number of target points in y direction for full channel calc

% Setting up channel

%            1 - top (perm)
%         +""""""""""""""""""""""""""""""+
%         |                              |
%  3-left |                              | 4-right
% (known) |                              | (known)
%         |                              |                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     |   
%         +""""""""""""""""""""""""""""""+
%         2 - bottom (perm)

% Channel Geometry
Lx = pi;
Ly = 1;
xmin = -Lx;
xmax = Lx;
ymin = -Ly;
ymax = Ly;
perm_factor=1;  % will make permeable section  -factor*Lx < x < factor*Lx
scalefactor = 0.05; % for the boundary velocity figures

% Discretization step
ds_x = (xmax - xmin)/(Ny1+1);
ds_y = (ymax - ymin)/(Ny2-1);

% Regularization parameter
ep =ds_y/20;%ds_y;
beta_value = Da*ep*0.04;

% Discretize Channel Walls (corners belong to left and right)
% top and bottom
stb = xmin+ds_x:ds_x:xmax-ds_x;
stb = stb';
% left and right
slr = ymin:ds_y:ymax;
slr = slr';

% Define coordinates on each wall (x-coord,y-coord)
% y1_top = stb;   y2_top = ymax*ones(size(stb)); % top wall coordinates (y1_top,y2_top)
% y1_bot = stb;   y2_bot = ymin*ones(size(stb)); % bottom wall coordinates (y1_bot,y2_bot)
% y1_left = xmin*ones(size(slr));   y2_left = slr; % left wall coordinates (y1_left,y2_left)
% y1_right = xmax*ones(size(slr));   y2_right = slr; % right wall coordinates (y1_right,y2_right)
% y1 = [y1_top; y1_bot; y1_left; y1_right]; %x-coordinates of all boundary points
% y2 = [y2_top; y2_bot; y2_left; y2_right]; %y-coordinates of all boundary points
% 
% wt = [ds_x*ones(size(y1_top)); ds_x*ones(size(y1_bot)); ds_y*ones(size(y1_left)); ds_y*ones(size(y1_right))];

% Define coordinates on each wall (x-coord,y-coord)
y1_top = stb;   y2_top = ymax*ones(size(stb)); % top wall coordinates (y1_top,y2_top)
y1_bot = flip(stb);   y2_bot = ymin*ones(size(stb)); % bottom wall coordinates (y1_bot,y2_bot)
y1_left = xmin*ones(size(slr));   y2_left = slr; % left wall coordinates (y1_left,y2_left)
y1_right = xmax*ones(size(slr));   y2_right = flip(slr); % right wall coordinates (y1_right,y2_right)
y1 = [y1_top; y1_right; y1_bot; y1_left]; %x-coordinates of all boundary points
y2 = [y2_top; y2_right; y2_bot; y2_left]; %y-coordinates of all boundary points
% 
wt = [ds_x*ones(size(y1_top)); ds_y*ones(size(y1_right)); ds_x*ones(size(y1_bot)); ds_y*ones(size(y1_left));];

% Exact velocities on boundaries using true solution (Bernardi, 2023)
 [u1_top_exact,u2_top_exact,p_top_exact] = permeablechannelexact(y1_top,y2_top,Da); %top velocity
 [u1_bot_exact,u2_bot_exact,p_bot_exact] = permeablechannelexact(y1_bot,y2_bot,Da); %bottom velocity
 [u1_left_exact,u2_left_exact,p_left_exact] = permeablechannelexact(y1_left,y2_left,Da); % left velocity
 [u1_right_exact,u2_right_exact,p_right_exact] = permeablechannelexact(y1_right,y2_right,Da); %right velocity

 u1_bd_exact = [u1_top_exact; u1_right_exact; u1_bot_exact; u1_left_exact]; %x-coordinates of all boundary velocities
 u2_bd_exact  = [u2_top_exact; u2_right_exact; u2_bot_exact; u2_left_exact]; %y-coordinates of all boundary velocities



% Indices of permeable region (where the boundary velocity is unknown)
%idx = find(-perm_factor*Lx<y1 & y1<perm_factor*Lx);
%beta = zeros(length(y1),1);
beta = zeros(Nseg,1);
idx = 1:(2*Ny1+2);
beta(idx) = beta_value;
% 
% % Create segment data matrix [x-idx, y-idx, ep, beta]
% % Order of segments: top (l->r), bottom (l->r), left (b->t), right (b->t)
% % segData = zeros(Nseg,4);
% % 
% % segData(:,1) = [2*Ny1+Ny2; (1:Ny1)'; 2*Ny1+1; (Ny1+1:2*Ny1)'; 
% %    (2*Ny1+1:2*Ny1+Ny2-1)'; (2*Ny1+Ny2+1:2*Ny1+2*Ny2-1)'];
% % segData(:,2) = [(1:Ny1)'; 2*Ny1+2*Ny2; (Ny1+1:2*Ny1)'; 2*Ny1+Ny2+1; 
% %    (2*Ny1+2:2*Ny1+Ny2)'; (2*Ny1+Ny2+2:2*Ny1+2*Ny2)'];
% % segData(:,3) = ep*ones(Nseg,1);
% % segData(:,4) = beta.*ones(Nseg,1);

% Going around channel instead
segData(:,1) = (1:Nseg)';
segData(:,2) = [(2:Nseg)'; 1];
segData(:,3) = ep*ones(Nseg,1);
segData(:,4) = beta.*ones(Nseg,1);

% Going around cc
segData(:,1) = (Nseg:-1:1)';
segData(:,2) = [(Nseg-1:-1:1)'; Nseg];

% Plotting segments
% figure
% hold on
% axis([-4,4,-1.5,1.5])
% for k = 1:Nseg
%     y = [y1(segData(k,1)) y2(segData(k,1));
%     y1(segData(k,2)) y2(segData(k,2))];
%     plot(y(:,1),y(:,2),'-','linewidth',3)
%     pause(0.2)
% end
% title('Segments')

% Define normal vectors
normals_top = zeros(length(stb),2); % unit normals for top:
normals_top(:,2) = 1;
normals_bot = zeros(length(stb),2); % unit normals for bottom
normals_bot(:,2) = -1;
normals_left = zeros(length(slr),2); % unit normals for left  wall 
normals_left(:,1) = -1;
normals_right = zeros(length(slr),2); % unit normals for right wall 
normals_right(:,1) = 1;
% normals on full boundary
normals = [normals_top; normals_right; normals_bot; normals_left];


% Points where velocity will be computed within the channel
xx1 = linspace(xmin+2*ds_x,xmax-2*ds_x,Nx1);
xx2 = linspace(ymin+2*ds_y,ymax-2*ds_y,Nx2);
dx_g = (xmax-xmin)/Nx1;
dy_g = (ymax - ymin)/Nx2;
[x1m, x2m] = ndgrid(xx1, xx2);
x1 = x1m(:); %x-coords of points where computing velocity
x2 = x2m(:); %y-coords of points where computing velocity
[x1gg, x2gg] = meshgrid(xx1, xx2);

%% Stokeslet and Source Doublet Velocities with Segments

I = eye(2,2);
x = [x1, x2]; % evaluation points
N = length(x1); % number of eval points
% Initialize segment velocities
u1seg_g = 0*x1; u2seg_g = 0*x2; 
u1sd_seg_g = 0*x1; u2sd_seg_g = 0*x2;
% Get forces from MRS code using exact boundary velocity
f = RegStokeslets2D_velocitytoforce([y1,y2], [y1,y2], [u1_bd_exact,u2_bd_exact], ep, mu,1,wt);
g = RegStokeslets2D_velocityto_gforce_permeable([y1,y2],[y1,y2],[u1_bd_exact,u2_bd_exact],ep,mu,1,idx, beta, normals, wt);
% Normals for segments
normalseg = normals(segData(:,1),:);

for m = 1:N % looping over evaluation points 
    for k = 1:Nseg % looping over segments
        % Current segment endpoints y(1,:) = yj(left), y(2,:) = yk(right)
        y = [y1(segData(k,1)) y2(segData(k,1));
            y1(segData(k,2)) y2(segData(k,2))];
        ep = segData(k,3);          
        [Q0,Q1,T02,T12,T22,T32] = st_segment_terms(x(m,:),y,ep);
        x0 = (x(m,:) - y(1,:))'; % x0 = x-yj
        v = ((y(1,:)-y(2,:)))'; % v = yj-yk (unnormalized tangent vector)
        L = norm(v); % segment length
        fj = f(segData(k,1),:);%/(L); % force at left endpoint of kth segment
        fk = f(segData(k,2),:);%/(L); % force at right endpoint of kth segment
        gj = g(segData(k,1),:);
        gk = g(segData(k,2),:);
        M2 = ((ep^2*T12 - Q1)*I + T12.*(x0*x0') + T22.*(x0*v' + v*x0') + T32.*(v*v'));
        M1 = ((ep^2*T02 - Q0)*I + T02.*(x0*x0')+ T12.*(x0*v' + v*x0') + T22.*(v*v')) - M2;
        u_temp = (M1*fj' + M2*fk')/(4*pi);
        u1seg_g(m) = u1seg_g(m)+ L*u_temp(1); 
        u2seg_g(m) = u2seg_g(m)+ L*u_temp(2);
        [T02,T04,T06,T12,T14,T16,T24,T26] = sd_segment_terms(x(m,:),y,ep);
        n = [normalseg(k,1); normalseg(k,2)];
        nx0 = n'*x0; 
        x0n = x0*n';
        vn = v*n';
        nn = n*n';
        D2 = ((T12 + ep^2*T14).*nn + nx0.*(-2*T14 -4*ep^2*T16).*x0n + nx0.*(- 2*T24 - 4*ep^2*T26).*vn);
        D1 = (T02 + ep^2*T04).*nn + nx0.*(-2*T04 - 4*ep^2*T06).*x0n + nx0.*(-2*T14 - 4*ep^2*T16).*vn  - D2;
        usd_temp = -(D1*gj' + D2*gk')/(2*pi); 
        u1sd_seg_g(m) = u1sd_seg_g(m) + L*usd_temp(1);
        u2sd_seg_g(m) = u2sd_seg_g(m) + L*usd_temp(2);
    end
end 

%% Now with function file
 % mu = 1;
 % x = [x1, x2]; % evaluation points
 % f = RegStokeslets2D_velocitytoforce([y1,y2], [y1,y2], [u1_bd_exact,u2_bd_exact], ep, mu,1,wt);
 % [Useg] = st_segments_forcetovelocity(x,[y1 y2],f,mu,segData);
%% Streamline Plot for Stokeslet Segment Velocity vs True Velocity

u1_seg_m = reshape(u1seg_g,size(xx1,2),size(xx2,2)); %x-coords of computed velocities
u2_seg_m = reshape(u2seg_g,size(xx1,2),size(xx2,2)); %y-coords of computed velocities
[uexact,vexact,~] = permeablechannelexact(x1gg,x2gg,Da);
u_error_seg = uexact - u1_seg_m';
v_error_seg = vexact - u2_seg_m';

figure
plot(y1, y2, 'k.')
hold on

% Computed streamlines
hh1_comp = streamline(x1gg, x2gg, u1_seg_m', u2_seg_m', x1gg(1:end, 1), x2gg(1:end, 1)); % streamlines starting at left wall
hh2_comp = streamline(x1gg, x2gg, -u1_seg_m', -u2_seg_m', x1gg(1:end, end), x2gg(1:end, end)); % streamlines starting at right wall
% Exact streamlines
hh1_exact = streamline(x1gg, x2gg, uexact, vexact, x1gg(1:end, 1), x2gg(1:end, 1)); % streamlines starting at left wall
hh2_exact = streamline(x1gg, x2gg, -uexact, -vexact, x1gg(1:end, end), x2gg(1:end, end)); % streamlines starting at right wall

strmLW = 2;
RGB = orderedcolors("gem");
set(hh1_comp, 'Color', RGB(1,:), 'linewidth', strmLW);
set(hh2_comp, 'Color', RGB(1,:), 'linewidth', strmLW);
set(hh1_exact, 'Color', RGB(3,:), 'LineStyle','--', 'linewidth', strmLW);
set(hh2_exact, 'Color', RGB(3,:), 'LineStyle','--', 'linewidth', strmLW);

title('Computed (blue) and exact (yellow) streamlines - Stokeslet Segments')
set(gca, 'FontSize', 14)
xlabel('$x$', 'Interpreter', 'latex')
ylabel('$y$', 'Interpreter', 'latex')
axis equal
axis([-4,4,-1.5,1.5])
%% MRS Stokeslet and Source Doublet Velocities for Comparisons

% NOTE: We are supposed to use midpoint rule here... maybe redo this with
% traprule?

% These aren't really comparable, especially for the source doublet
% calculation since we don't have a true solution for comparison.

 ep = ds_y;

% Get forces from MRS code using exact boundary velocity
f = RegStokeslets2D_velocitytoforce([y1,y2], [y1,y2], [u1_bd_exact,u2_bd_exact], ep, mu,1,wt);
%g = RegStokeslets2D_velocityto_gforce_permeable([y1,y2],[y1,y2],[u1_bd_exact,u2_bd_exact],ep,mu,1,idx, beta, normals, wt);

ep = ds_y;
[ug_st_check] = RegStokeslets2D_forcetovelocity([y1,y2],f,[x1,x2],ep,mu,1,wt);
ug_st = ug_st_check(:,1);
vg_st = ug_st_check(:,2);
%u_seg_g = [u1seg_g u2seg_g];

[ug_sd_check] = RegStokeslets2D_gtovelocity([y1,y2],g,[x1,x2],ep,mu,1,beta, normals, wt);
 ug_sd = ug_sd_check(:,1);
 vg_sd = ug_sd_check(:,2);
usd_seg_g = [u1seg_g u2seg_g];
%% Streamline Plot for MRS Stokeslet Velocity vs True Velocity

um_st = reshape(ug_st,size(xx1,2),size(xx2,2)); %x-coords of computed velocities
vm_st = reshape(vg_st,size(xx1,2),size(xx2,2)); %y-coords of computed velocities

u_error_mrs = uexact - um_st';
v_error_mrs = vexact - vm_st';

% Plot streamlines comparing exact vs computed velocities using Stokeslets
figure
plot(y1, y2, 'k.')
hold on

% Computed streamlines
hh1_comp = streamline(x1gg, x2gg, um_st', vm_st', x1gg(1:end, 1), x2gg(1:end, 1)); % streamlines starting at left wall
hh2_comp = streamline(x1gg, x2gg, -um_st', -vm_st', x1gg(1:end, end), x2gg(1:end, end)); % streamlines starting at right wall
% Exact streamlines
hh1_exact = streamline(x1gg, x2gg, uexact, vexact, x1gg(1:end, 1), x2gg(1:end, 1)); % streamlines starting at left wall
hh2_exact = streamline(x1gg, x2gg, -uexact, -vexact, x1gg(1:end, end), x2gg(1:end, end)); % streamlines starting at right wall

strmLW = 2;
RGB = orderedcolors("gem");
set(hh1_comp, 'Color', RGB(1,:), 'linewidth', strmLW);
set(hh2_comp, 'Color', RGB(1,:), 'linewidth', strmLW);
set(hh1_exact, 'Color', RGB(3,:), 'LineStyle','--', 'linewidth', strmLW);
set(hh2_exact, 'Color', RGB(3,:), 'LineStyle','--', 'linewidth', strmLW);

title('Computed (blue) and exact (yellow) streamlines - MRS')
set(gca, 'FontSize', 14)
xlabel('$x$', 'Interpreter', 'latex')
ylabel('$y$', 'Interpreter', 'latex')
axis equal
axis([-4,4,-1.5,1.5])

% MRS Errors
[MRS_Stokeslet_errors] = calculate_errors(u_error_mrs,v_error_mrs,xx1,xx2,dx_g,dy_g, 'MRS Stokeslet Only')


% Segment Errors
[Segment_Stokeslet_errors] = calculate_errors(u_error_seg,v_error_seg,xx1,xx2,dx_g,dy_g, 'Segment Stokeslet Only')


%% Comparing Source Doublet Segment Velocity to MRS Source Doublet Velocity
u_sd_seg_m = reshape(u1sd_seg_g,size(xx1,2),size(xx2,2)); %x-coords of computed velocities
v_sd_seg_m = reshape(u2sd_seg_g,size(xx1,2),size(xx2,2)); %y-coords of computed velocities

u_sd_st_m = reshape(ug_sd,size(xx1,2),size(xx2,2)); %x-coords of computed velocities
v_sd_st_m = reshape(vg_sd,size(xx1,2),size(xx2,2)); %y-coords of computed velocities

usd_diff_mrs = u_sd_seg_m' - u_sd_st_m';
vsd_diff_mrs = v_sd_seg_m' - v_sd_st_m';

[MRS_Stokeslet_errors] = calculate_errors(usd_diff_mrs,vsd_diff_mrs,xx1,xx2,dx_g,dy_g, '(Segment-MRS) differences SD Velocity')

%% Exact solutions
function [ug,vg,pg] = permeablechannelexact(x,y,Da)

% Newton's method to find lambda1
L1=sqrt(2*Da);
f=@(L) tan(L)^2-tan(L)/L+1-2*Da;
L1 = fzero(f,L1);
%disp(['lambda_1 = ',num2str(L1)])

% Compute g and g' from FB paper
g  = cot(L1)*(Da-0.5)   *sin(L1*y)+0.5*y.*cos(L1*y);
gp = cot(L1)*(Da-0.5)*L1*cos(L1*y)+0.5*cos(L1*y)-0.5*L1*y.*sin(L1*y);

% Compute exact solution using g and g'
C = -4;
pg = C*sinh(L1*x) .* cos(L1*y);
ug = -C/L1 * cosh(L1*x) .* gp;
vg =  C    * sinh(L1*x) .* g ;

end
%% Errors
function [errors] = calculate_errors(u_error,v_error,xx1,xx2,dx,dy, method_name)

u_error_max = max(max(abs(u_error)));
v_error_max = max(max(abs(v_error)));
u_error_L1 = dx*dy*sum(sum(abs(u_error)));
v_error_L1 = dx*dy*sum(sum(abs(v_error)));
u_error_L2 = sqrt(dx*dy*sum(sum(u_error.^2)));
v_error_L2 = sqrt(dx*dy*sum(sum(v_error.^2)));

ColNames = {'Max','L1','L2'};
Emax = [u_error_max; v_error_max];
EL1 = [u_error_L1; v_error_L1];
EL2 = [u_error_L2; v_error_L2];
errors = table(Emax, EL1, EL2,'VariableNames',ColNames);

figure
subplot(2,1,1);
surf(xx1, xx2, u_error, 'EdgeColor', 'none');
colorbar;
title(['Error in u (horizontal) velocity -- ' method_name]);
xlabel('x');
ylabel('y');
axis equal;
view(2); % View from above for 2D representation

subplot(2,1,2);
surf(xx1, xx2, v_error, 'EdgeColor', 'none');
colorbar;
title(['Error in v (vertical) velocity -- ' method_name]);
xlabel('x');
ylabel('y');
axis equal;
view(2); % View from above for 2D representation

end