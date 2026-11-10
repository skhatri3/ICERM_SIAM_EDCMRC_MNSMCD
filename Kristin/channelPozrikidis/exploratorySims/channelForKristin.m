
clear all
 figure(1)
%% Parameters to set

% Model Parameters
Da = 0.4; % Darcy Number
Dab = Da; % Darcy Number
mu = 1; % Viscosity
% Channel Geometry -(Lx+c0) < x < Lx+c0, -H < h < H
H=1; Ly = H;
Lx = 1.5;  c0 = 2;
% The inlet flow is parabolic based on a constant pressure gradient -G
G=4;   ep = H/14/2;
% Find lambda by Newton's method
F=@(z) z*(Da/H-1)+Da/H*z.*cos(z)+sin(z);
Fp=@(z) (Da/H-1)+Da/H*cos(z)-Da/H*z*sin(z)+cos(z);
z=1; for iter=1:16
z=z-F(z)/Fp(z); end 
lam = z/2/Ly;
%Pozrikidis version uses alfa and gam
alfa = sqrt( 3*mu*Dab/H^3 );
gam = -G/alfa/cosh(alfa*Lx);
% define the functions g(y) and g'(y) from the notes
g=@(y) (Da - H/2)*cot(lam*H)*sin(lam*y) + y/2.*cos(lam*y);
gp=@(y) (Da - H/2)*lam*cot(lam*H)*cos(lam*y) + 1/2.*cos(lam*y) - lam*y/2.*sin(lam*y);
% compute the coeficients/parameters
%slopem = (alfa*gam*cosh(lam*Lx)+G)/(1-cosh(lam*Lx));
slopem = 0;
% slopem = 0;
bR = Lx*(slopem-G) *(tanh(lam*Lx)/lam/Lx -1 );
c1=(slopem-G)/lam/cosh(lam*Lx)/cos(lam*H) ; 
bL = - bR;

% Francesca's version
         % p = c1 * cos(lam*y) * sinh(lam*x) - slopem*x
         % slopem is included to match Pozrikidis

 % for | x | > Lx (solid boundaries), p = -G*x+bL and p = -G*x+bR
% bR = Lx*(slopem-G) *(tanh(lam*Lx)/lam/Lx -1 );
% c1=(slopem-G)/lam/cosh(lam*Lx)/cos(lam*H) ; 
% bL = - bR;

xmin = -Lx-c0;
xmax = Lx+c0;
ymin = -H;
ymax = H;

influx = G*H^3*2/3;
topflux = Da*(G-slopem)*(cosh(lam*Lx)-1)/lam^2/cosh(lam*Lx) + Da*slopem*Lx^2/2;
disp(['Inlet flux = ',num2str(influx)])
disp(['Top flux in half-channel = ',num2str( topflux)])
disp(['Bottom flux in half-channel = ',num2str( topflux)])
disp(['Therefore, the flux out of the channel at the midcross-section is ',num2str(influx-2*topflux)])

xx = (-Lx-c0 : 2*(Lx+c0)/100 : Lx+c0)';
yy = (-H : H/14 : H)';
%pressures at (x,y)=(Lx,H)
p1L = c1 * sinh(lam*Lx) *cos(lam*H) - Lx*slopem; %Francesca
vtop = c1 * sinh(lam*xx) *g(H) - Da*xx*slopem; %Francesca
p2L = gam*sinh(alfa*Lx);                                           % Pozrikidis
disp('pressures'),disp([p1L p2L])
disp('lambda, alpha c1 gam'),disp([lam alfa c1 gam])
P1 = (c1 * sinh(lam*xx) *cos(lam*H) - xx*slopem) .*(xx<Lx).*(xx>-Lx) ...
    +(-G*(xx-Lx)+p1L).*(xx>=Lx)  +(-G*(xx+Lx)-p1L).*(xx<=-Lx);
P2 = gam*sinh(alfa*xx).*(xx<Lx).*(xx>-Lx) ...
    +(-G*(xx-Lx)+p2L).*(xx>=Lx)  +(-G*(xx+Lx)-p2L).*(xx<=-Lx); 

% All four sides of channel
xx2 = xx(2:end-1);
y1 = [xx2 ; xx2 ; -(Lx+c0)*ones(size(yy)) ;  (Lx+c0)*ones(size(yy))];
y2 = [H*ones(size(xx2)) ; -H*ones(size(xx2)) ; yy ; yy];
u1 = [0*ones(size(xx2)) ; 0*ones(size(xx2)) ; G/2*(H^2-yy.^2)  ; G/2*(H^2-yy.^2)];
vtop2 = 0*xx2 + (c1 * sinh(lam*xx2) *g(H) - Da*xx2*slopem).*(abs(xx2)<Lx); %Francesca
u2 = [vtop2 ; -vtop2 ; 0*ones(size(yy)) ; 0*ones(size(yy))];
%ftemp = RegStokeslets2D_velocitytoforce([y1,y2], [y1,y2], [u1,u2], ep, mu, 1,ones(size(y1)));
%----
Uall = [u1'; u2']; Uall = Uall(:);
StokesA = MatrixStokeslet2D([y1,y2], [y1,y2],ep,1);
ftmpall = StokesA\Uall; ftemp = [ftmpall(1:2:end) ftmpall(2:2:end)];
[xg,yg]=meshgrid( xx , (-H : H/14 : H)); nx = length( (-H : H/14 : H) ); ny = length(xx);
[ug,vg] = Stokeslet2Dftou([xg(:),yg(:)],[y1(:),y2(:)],ftemp,ep,1);
ug = reshape(ug,size(xg));
vg = reshape(vg,size(xg));

% Compute pressure from forces
addpath('../stokeslets_codes'); 
blob_num = 2;
ds_x = xx2(2) - xx2(1);
ds_y = yy(2 ) -yy(1);
wt = [ds_x*ones(size(xx2)); ds_x*ones(size(xx2)); ds_y*ones(size(yy)); ds_y*ones(size(yy))];
P = RegStokeslets2D_forcetopressure([y1,y2], ftemp, [xg(:), yg(:)], ep, mu, blob_num, wt);
pg_st = reshape(P, size(xg));
% plot pressure
figure;
pcolor(xg, yg, pg_st);
shading interp;
axis equal
colorbar('southoutside', 'FontSize', 16)
ax = gca; ax.FontSize = 16;

%%
idx=1:1:101;
figure(2),clf, plot(y1,y2,'k.')%quiver(xg,yg,ug,vg,1);
axis equal, hold on
streamline(xg(:,idx) ,yg(:,idx),-ug(:,idx),-vg(:,idx),xg(:,51),yg(:,51));  hold on
streamline(xg(:,idx),yg(:,idx),ug(:,idx),vg(:,idx),xg(:,51),yg(:,51));
h=streamline(xg(:,idx),yg(:,idx),-ug(:,idx),-vg(:,idx),xg(:,end),yg(:,end)); set(h,'Color','Red')
h=streamline(xg(:,idx),yg(:,idx),ug(:,idx),vg(:,idx),xg(:,1),yg(:,1));  set(h,'Color','Red')
pause(2),set(gca,'FontSize',16),axis([-Lx-c0-0.5 Lx+c0+0.5 -1.5 1.5])
quiver(y1,y2,u1/10,u2/10,0,'b')
xlabel(['Da = ', num2str(Da),  ', ', ...
    'G = ', num2str(G),  ', ', ...
    '(L,H) = (', num2str(Lx),',',num2str(H),')',  ', ', ...
    'c = ', num2str(c0),   ...
    ''])

%
xxx = [(1e-6 : 1e-7 : 1e-4) (2e-4 : 1e-4 : 1e-2) (2e-2 : 1e-2 : 0.3 )]';
figure(1),hold off
plot(xx,P1-xx*0,'b-' ,xx,P2 ,'r--','LineWidth',2)
grid on, title(['Pressure at the top boundary, m = ',num2str(slopem)])
legend('Francesca','Pozrikidis')
set(gca,'FontSize',16)
[lam alfa gam G Da]
%%
figure(3)
pg=0*xg; pg(1,:)=P1; pg(end,:)=P1; pg(:,1)=P1(1); pg(:,end)=P1(end);
pg = ones(nx,1)*P1';
%levls = [(12:-1:3) 2.43 2.17 1.91 1.66 1.4 1.16 0.9 0.63 0.38 0.12]; levls=[levls -levls];
levls = [(12:-1:0) ]; levls=[levls -levls];
dx = mean(diff(xx)); dy = mean(diff(yy));
for iter=1:4001
    tmp = (pg(1:nx-2,2:ny-1) + pg(3:nx,2:ny-1))*dx^2  + (pg(2:nx-1,1:ny-2) + pg(2:nx-1,3:ny) )*dy^2;
    tmp = tmp/2/(  dx^2 + dy^2 );
    if mod(iter,600)==1,
            sum(sum( abs(tmp-pg(2:nx-1,2:ny-1) ).^2  *dx*dy )    );
            contourf(xg,yg,pg, levls,'k'),axis equal,pause(0.2)
    end
       corr= sqrt(sum(sum( abs(tmp-pg(2:nx-1,2:ny-1) ).^2  *dx*dy )   ) );
        pg(2:nx-1,2:ny-1) = tmp;
end
disp(['the L2 norm of the last correction was ',num2str(corr)])
% levls = [(12:-1:3) 2.43 2.17 1.91 1.66 1.4 1.16 0.9 0.63 0.38 0.12]; levls=[levls -levls];
[C,h] =contourf(xg,yg,pg, levls);axis equal
clabel(C,h)
title('Contours of Francesca''s pressure')
xlabel(['Da = ', num2str(Da),  ', ', ...
    'G = ', num2str(G),  ', ', ...
    '(L,H) = (', num2str(Lx),',',num2str(H),')',  ', ', ...
    'c = ', num2str(c0),   ...
    ''])
set(gca,'FontSize',16)
return





function A = MatrixStokeslet2D(x,y,d,Mu)

%
% x = evaluation points (matrix of size Nx3)
% y = location of the Stokeslets (matrix of size Nx3)
% d = delta
%

d2 = d^2;

nx = length(x(:,1));
ny = length(y(:,1));


% So x is nx by 3, so is y

A = zeros(2*nx,2*ny);

for i = 1 : ny

      columnid = 2*(i-1);
      columnid1 = columnid+1;
      columnid2 = columnid+2;

      dx1 = x(:,1)-y(i,1);
      dx2 = x(:,2)-y(i,2);
      
      r2  = dx1.^2 + dx2.^2 + d2;      
      
      H1 = -0.5*log(r2) + d2./r2;
      H2 = 1./r2;
      

      A(1:2:end,columnid1) = H1 + dx1.*dx1.*H2;
      A(1:2:end,columnid2) = dx1.*dx2.*H2 ;

      A(2:2:end,columnid1) = dx2.*dx1.*H2 ;
      A(2:2:end,columnid2) = H1 + dx2.*dx2.*H2 ;

end
A = A/(4*pi*Mu);

end%MatrixStokeslet2D
%%%


function [u1,u2] = Stokeslet2Dftou(x,y,f,d,Mu)

%
% x = evaluation points (matrix of size Nx3)
% y = location of the Stokeslets (matrix of size Nx3)
% d = delta
%

d2 = d^2;

nx = length(x(:,1));
ny = length(y(:,1));


% So x is nx by 2, so is y

u1 = zeros(nx,1);
u2 = zeros(nx,1);

for k = 1 : ny

  
      dx1 = x(:,1)-y(k,1);
      dx2 = x(:,2)-y(k,2);
      
      r2  = dx1.^2 + dx2.^2 + d2;      
      H1 = -0.5*log(r2) + d2./r2;
      H2 = 1./r2;
      
     fdotx = f(k,1)*dx1 + f(k,2)*dx2;
      u1= u1 + H1*f(k,1) + fdotx.*dx1.*H2;
      u2= u2 + H1*f(k,2) + fdotx.*dx2.*H2;

end
u1 = u1/(4*pi*Mu);
u2 = u2/(4*pi*Mu);

end%MatrixStokeslet2D
%%%