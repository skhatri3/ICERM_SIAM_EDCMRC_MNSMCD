%Example 3 of Cortez, Fluids 2021
%Channel with inflow and part of membrane permeable

%Developed by Ricardo Cortez, Brittany Leathers, and Michaela Kubacki
%July 2024

%Updated September 2024 to allow for different choices of epsilon for
%Stokeslets and Source Doublets.

% Looking at channel without permeability.

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
Nvals = [10 20 40 80 160 320];


% % C1 value range
% C1vals = linspace(0.1,4,20);

% % Keeping ep2 constant (looking for C1 in ep1 = C1ds)
% %ep2 = 0.0002;
% ep2 = 0.05/sqrt(5);

%blob choice (1 is phi and 2 is psi)
blob=2;

%constant determining inflow velocity profile
a=4;
 
%channel set-up
Lx = 5;
Ly = 1;
xmin = 0;
xmax = xmin + Lx;
ymin = 0;
ymax = ymin + Ly;
perm_min=5/3;  %start of permeable part
perm_max=10/3;  %end of permeable part

% Set up error matrix and initialize figure
%ErrorMatrix = zeros(length(Nvals),length(epvals));
% Net Flow Error Plot
fig1 = figure;
hold on;
set(gca, 'FontSize', 17);
xlabel('$\epsilon$', 'FontSize', 18,'Interpreter','latex')
ylabel('Error', 'FontSize', 17)
title('Horizontal Velocity Error (L2), non-permeable, using blob $\psi$',...
    'FontSize', 18,'Interpreter','latex')
hold off;


fig2 = figure;
hold on;
set(gca, 'FontSize', 17);
xlabel('$\epsilon$', 'FontSize', 18,'Interpreter','latex')
ylabel('Error', 'FontSize', 17)
title('Horizontal Velocity Error (Linf), non-permeable, using blob $\psi$',...
    'FontSize', 18,'Interpreter','latex') 
hold off;

for k = 1:length(Nvals)
    N = Nvals(k);
    %discretization of channel
    ds = (ymax-ymin)/N;
    epvals = logspace(-4,0,20);
    % usolutions=cell(1,length(N)); %for storing u solutions to compare
    % vsolutions=cell(1,length(N));

    for i = 1:length(epvals)
        ep1 = epvals(i);

        %permeability coefficient (assuming constant for the permeable region) 
        % for blob 1
        %b = 0.011843080738230*ep2; % for b = 0.0001 (beta 1)
        %b = 0.004247164369662*ep2; % for b = 0.00018
        %b = 0.0282306994128518*ep2; % for b = 0.00024 flow


        %permeability coefficient (assuming constant for the permeable region) 
        % for blob 2
        % b=0.010733*ep2; % for b = 0.00024 (beta 3)
        % b = 0.00805*ep2; % for b = 0.00018 (beta 2)
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

        y2mid = y2right;
        y1mid = xmax/2*ones(size(y2_side));

        y1f=[y1; y1right];
        y2f=[y2; y2right];
        %indices of permeable region
        I=find(y1_top<2/3*Lx & y1_top>1/3*Lx);
        I2=find(y1_top>=2/3*Lx | y1_top<=1/3*Lx);
        %unit normals
        normals=[normals_top; normals_bot; normals_side];
        %beta vector for function;
        beta= zeros(length(y1),1);

 % TURN OFF PERMEABILITY (beta = 0)
 %       beta(I)=b; 

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

% IF TURN OFF PERMEABILITY, NO LONGER NECESSARY, since this is just a Stokeslet problem
%        %compute g
%        g=RegStokeslets2D_velocityto_gforce_permeable_diff_ep([y1,y2],[y1,y2],...
%        [u1,u2],ep1,ep2,mu, blob, I, beta, normals);

%        %Find velocity in permeable region:
%        y1b=y1(I);
%        y2b=y2(I);

%        [u_beta] = RegStokeslets2D_gtovelocity([y1,y2],g, [y1b,y2b],ep2,mu,blob, beta, normals);

%        %Put in the new velocities for the permeable part
%        %permeable part : temporary velocity for finding g
%        u1(I)=u_beta(:,1);
%        u2(I)=u_beta(:,2);

        %computing the force
        f = RegStokeslets2D_velocitytoforce([y1,y2],[y1,y2],[u1,u2],ep1,mu, blob);
        f1 = f(:,1);
        f2 = f(:,2);  

        % Calculate velocities on right hand side too
        Ufull=RegStokeslets2D_forcetovelocity([y1,y2],[f1,f2],[y1f,y2f],ep1,mu, blob, ds);
        ufull=Ufull(:,1); vfull=Ufull(:,2);

        % Calculus velocities in middle
        Umid=RegStokeslets2D_forcetovelocity([y1,y2],[f1,f2],[y1mid,y2mid],ep1,mu, blob, ds);
        umid=Umid(:,1); vmid=Umid(:,2);


        %calculate on grid
        [xgg,ygg] = meshgrid(ds:4*ds:Lx-ds, ds:4*ds:Ly-ds); xg=xgg(:); yg=ygg(:);

        Ugrid=RegStokeslets2D_forcetovelocity([y1,y2],[f1,f2],[xg,yg],ep1,mu, blob, ds);
        ug=reshape(Ugrid(:,1), length(ds:4*ds:Ly-ds), length(ds:4*ds:Lx-ds)); 
        vg=reshape(Ugrid(:,2), length(ds:4*ds:Ly-ds), length(ds:4*ds:Lx-ds)); 
        speed=sqrt(ug.^2+vg.^2);

        % Plot figure
        % sk=2;
        % figure;%subplot(211)
        % plot(y1f,y2f,'k.'),hold on
        % % quiver(xgg,ygg,ug,vg,0.8,'r','LineWidth',1)
        % quiver(y1f(1:6*sk:end),y2f(1:6*sk:end),ufull(1:6*sk:end),vfull(1:6*sk:end),0,'k','LineWidth',2)
        % % quiver(xe(1:8:end),ye(1:8:end),usuck(1:8:end),vsuck(1:8:end),0,'r')
        % surf(xgg,ygg,-speed),view(2),shading interp
        % hh1=streamline(xgg,ygg,ug ,vg ,xgg(1:end,1 ),ygg( 1:end,1 ));
        % hh2=streamline(xgg,ygg,ug ,vg ,xgg(1:end,end),ygg( 1:end,end));
        % set(hh1,'Color','black');hold off,axis equal,
        % set(hh2,'Color','black');hold off,axis equal,axis([-0.10 5.5 -0.55 1.75  ])
        % title(['\epsilon= ',num2str(ep1), ', N = ',num2str(N)])
        % % clim([0 1]);
        % colorbar('Ticks',[-1 , -0.8, -0.6,-0.4,-0.2,0 ],...
        %     'TickLabels',{'1','0.8','0.6','0.4','0.2','0'},...
        %     'Direction','reverse')
        % hold off  
        % set(gca, 'FontSize', 16)
        % fig_title = ['epsilon= ',num2str(ep1), 'N=',num2str(N), '.fig'];
        % savefig([currPath,'/',fig_title])

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

        % Change to evaluate at x = 2.5 instead of at outlet
        %x2_out= (ymin+ds/2:ds:ymax-ds/2)';
        %x1_out= xmax/2*ones(size(x2_out)); % xmax/2 for x = 2.5
        %uu = RegStokeslets2D_forcetovelocity([y1,y2],[f1,f2],...
        %    [x1_out,x2_out],ep1,mu, blob);
        %Rout=-ds*sum(dot(normals_side, uu));

        % Calculus velocities in middle
        Umid=RegStokeslets2D_forcetovelocity([y1,y2],[f1,f2],[y1mid,y2mid],ep1,mu, blob, ds);
        umid=Umid(:,1); vmid=Umid(:,2);

        diffu = Umid(:,1)-u1_side;
        diffv = Umid(:,2)-u2_side;
        uL2error(i) = sqrt(ds^2*sum(diffu.^2));
        vL2error(i) = sqrt(ds^2*sum(diffv.^2));
        uLinferror(i) = max(max(abs(diffu)));
        vLinferror(i) = max(max(abs(diffv)));
        %error(i) = norm()
        %error(i)=Rin+Rtop+Rout;

    end


    % Save Current Error results and update figure
    %ErrorMatrix(k,:) = abs(error);

    name = ['N = ' num2str(N)];
    
    figure(fig1); hold on;
    loglog(epvals, uL2error,'o-', 'LineWidth', 2.5, 'MarkerSize', ...
       10,'DisplayName',name)
    hold off;

    figure(fig2); hold on;
    loglog(epvals,uLinferror,'x--','Linewidth',2.5, 'MarkerSize',10,'DisplayName',name)
    hold off;

 
end

figure(fig1)
set(gca,'YScale','log');
set(gca,'xscale','log');
legend


fig_title = ['Horizontal Error (L2) - no permeability','.fig'];
savefig([currPath,'/',fig_title])

figure(fig2)
set(gca,'YScale','log');
set(gca,'xscale','log');
legend



fig_title = ['Horizontal Error (Linf)- no permeability','.fig'];
savefig([currPath,'/',fig_title])
