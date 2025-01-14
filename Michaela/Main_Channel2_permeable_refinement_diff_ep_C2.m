%Example 3 of Cortez, Fluids 2021
%Channel with inflow and part of membrane permeable

%Developed by Ricardo Cortez, Brittany Leathers, and Michaela Kubacki
%July 2024

%Updated September 2024 to allow for different choices of epsilon for
%Stokeslets and Source Doublets. 

clear all 
close all

currDate = datestr(datetime);

mkdir(currDate)
currPath = [pwd,'/',currDate];

format long

%% Parameters to set 

% setting the viscosity
mu = 1; 

% number of points on boundary where velocity is set and force is computed 
Nvals = [10 20 40 80 160];

% C2 values
%C2vals = 0.01:0.005:0.1;
C2vals = logspace(-3,0,20);

% Keep ep1 constant (looking for C2 in ep2 = C2 (ds)^(1/3))
ep1 = 0.0002;

%blob choice
blob=2;

%constant determining inflow velocity profile
a=4;
 
%channel
Lx = 5;
Ly = 1;
xmin = 0;
xmax = xmin + Lx;
ymin = 0;
ymax = ymin + Ly;
perm_min=5/3;  %start of permeable part
perm_max=10/3;  %end of permeable part

% Set up error matrix and initialize figure
ErrorMatrix = zeros(length(Nvals),length(C2vals));
% Net Flow Error Plot
figure;
hold on;
set(gca, 'FontSize', 17);
xlabel('$C_2$', 'FontSize', 18,'Interpreter','latex')
ylabel('Error', 'FontSize', 17)
title('Net Flow Rate, $\beta = 0.00018$, using blob $\psi$',...
    'FontSize', 18,'Interpreter','latex')  
for k = 1:length(Nvals)
    N = Nvals(k); 

    % usolutions=cell(1,length(N)); %for storing u solutions to compare
    % vsolutions=cell(1,length(N));

    for i = 1:length(C2vals)
        C2=C2vals(i);
   
        %discretization of channel
        ds = (ymax-ymin)/N;

        %regularization parameter
        ep2 = C2 * ds^(1/3);

        %permeability coefficient (assuming constant for the permeable region) 
        % for blob 1
        %b = 0.011843080738230*ep2; % for b = 0.0001 (beta 1)
        %b = 0.004247164369662*ep2; % for b = 0.00018
        %b = 0.0282306994128518*ep2; % for b = 0.00024 flow
        
        
        %permeability coefficient (assuming constant for the permeable region) 
        % for blob 2
        % b=0.010733*ep2; % for b = 0.00024 (beta 3)
         b = 0.00805*ep2; % for b = 0.00018 (beta 2)
        % b = 0.004472*ep2; % for b = 0.0001 (beta 1)

        %discretization of top and bottom:
        s = ds/2:ds:xmax-ds/2;
        s = s';
        %top wall coordinates (x,y) = (y1_top,y2_top)
        y1_top=s;
        y2_top = ymax*ones(size(y1_top));
        %unit normals for top: 
        normals_top=zeros(length(y1_top),2);
        normals_top(:,2)=1;
        %bottom wall coordinates (x,y) = (y1_bot,y2_bot)
        y1_bot = s;
        y2_bot = ymin*ones(size(y1_bot));
        %unit normals for bottom
        normals_bot=zeros(length(y1_bot),2); 
        normals_bot(:,2)=-1;
        %left-hand wall coordinates (x,y) = (y1_side,y2_side)
        y2_side = (ymin+ds/2:ds:ymax-ds/2)';
        y1_side = xmin*ones(size(y2_side));
        %unit normals for side
        normals_side=zeros(length(y1_side),2);
        normals_side(:,1)=-1;
        
        %coordinates of boundary points
        y1 = [y1_top;y1_bot; y1_side];
        y2 = [y2_top;y2_bot; y2_side];
        %For use at very end:
        y2right = (ymin+ds/2:ds:ymax-ds/2)';
        y1right = xmax*ones(size(y2_side));
        y1f=[y1; y1right];
        y2f=[y2; y2right];
        %indices of permeable region
        I=find(y1_top<2/3*Lx & y1_top>1/3*Lx);
        I2=find(y1_top>=2/3*Lx | y1_top<=1/3*Lx);
        %unit normals
        normals=[normals_top; normals_bot; normals_side];
        %beta vector for function;
        beta= zeros(length(y1),1); 
        beta(I)=b;

        %velocity on boundary
        %top
        u1_top=zeros(size(y1_top));
        u2_top=zeros(size(y1_top));
        %bottom
        u1_bot=zeros(size(y1_bot));
        u2_bot=zeros(size(y1_bot));
        %side: poiseulle flow
        u1_side=a*(y2_side/Ly.*(1-y2_side/Ly));
        u2_side=zeros(size(y1_side));
        
        u1 = [u1_top; u1_bot; u1_side];
        u2 = [u2_top; u2_bot; u2_side];
         
        %
        %compute g
        g=RegStokeslets2D_velocityto_gforce_permeable_diff_ep([y1,y2],[y1,y2],...
    [u1,u2],ep1,ep2,mu, blob, I, beta, normals);

        %Find velocity in permeable region:
        y1b=y1(I);
        y2b=y2(I);

        [u_beta] = RegStokeslets2D_gtovelocity([y1,y2],g, [y1b,y2b],ep2,mu,blob, beta, normals);

        
        %Put in the new velocities for the permeable part
        %permeable part : temporary velocity for finding g
        u1(I)=u_beta(:,1);
        u2(I)=u_beta(:,2);

        %computing the force 
        f = RegStokeslets2D_velocitytoforce([y1,y2],[y1,y2],[u1,u2],ep1,mu, blob);
        f1 = f(:,1);
        f2 = f(:,2);  
        
        % Calculate velocities on right hand side too
        Ufull=RegStokeslets2D_forcetovelocity([y1,y2],[f1,f2],[y1f,y2f],ep1,mu, blob, ds);
        ufull=Ufull(:,1); vfull=Ufull(:,2);
        
        %calculate on grid
        [xgg,ygg] = meshgrid(ds:4*ds:Lx-ds, ds:4*ds:Ly-ds); xg=xgg(:); yg=ygg(:);
        Ugrid=RegStokeslets2D_forcetovelocity([y1,y2],[f1,f2],[xg,yg],ep1,mu, blob, ds);
        ug=reshape(Ugrid(:,1), length(ds:4*ds:Ly-ds), length(ds:4*ds:Lx-ds)); 
        vg=reshape(Ugrid(:,2), length(ds:4*ds:Ly-ds), length(ds:4*ds:Lx-ds)); 
        speed=sqrt(ug.^2+vg.^2);

        %calculate on grid
        % dx=1/160;
        % dx_s(i)=dx;
        % Nx=round(xmax/dx);
        % Ny=round(ymax/dx);
        % xg=dx*(0:Nx-1)+xmin;
        % yg=dx*(0:Ny-1)+ymin;
        % [xg,yg]=ndgrid(xg,yg);
        % xgv=reshape(xg, Nx*Ny,1);
        % ygv=reshape(yg, Nx*Ny,1);
        % Ugrid=RegStokeslets2D_forcetovelocity([y1,y2],[f1,f2],[xgv,ygv],ep1,mu, blob, ds);
        % ug=reshape(Ugrid(:,1), Nx, Ny); 
        % vg=reshape(Ugrid(:,2), Nx, Ny); 
        % speed=sqrt(ug.^2+vg.^2);
        % usolutions{i}=ug;
        % vsolutions{i}=vg;
        % 
        % xgo=xg';
        % ygo=yg';
        
        % Check flow rates 
        
        % inlet
        Rin=ds*sum(dot(normals_side, [u1_side u2_side]));
        
        %top
        normals_top=zeros(length( y1_top),2);
        normals_top(:,2)=1;
        Rtop=ds*sum(dot(normals_top, [ufull(1:length(y1_top)), vfull(1:length(y1_top))]));
        
        %outlet
        x2_out= (ymin+ds/2:ds:ymax-ds/2)';
        x1_out= xmax*ones(size(x2_out));
        uu = RegStokeslets2D_forcetovelocity([y1,y2],[f1,f2],...
            [x1_out,x2_out],ep1,mu, blob);
        Rout=-ds*sum(dot(normals_side, uu));
        
        error(i)=Rin+Rtop+Rout;

    end

    % Save Current Error results and update figure
    ErrorMatrix(k,:) = abs(error);   
    name = ['N = ' num2str(N)];
    loglog(C2vals, ErrorMatrix(k,:),'o-', 'LineWidth', 2.5, 'MarkerSize', ...
           10,'DisplayName',name)
     
end

hold off
set(gca,'XScale','log','YScale','log');
legend


fig_title = ['Finding C2','.fig'];
savefig([currPath,'/',fig_title])