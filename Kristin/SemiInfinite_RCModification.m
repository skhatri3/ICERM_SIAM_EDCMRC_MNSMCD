%function SemiInfinite_RCModification
clear
figure(2)
  NNh = [  60 ] ; %do various runs with NNh points per inlet unit length
  hh0 = 1./NNh;
  dd = hh0*1.0; % blob size
  netflx=0*ones(size(NNh));
jj=1

Da = 0.3; % Darcy number

% for jj=1:length(hh0)
%     jj;
% flag = 0; % 0 = do not use jumps;   1 = use jumps
% if flag==1, clr = 'r.-'; else clr = 'b.-'; end


% Set up the channel  -L/2<x<L/2, -H/2<y<H/2
% target discretization size h0
% 
%                    1
%    +""""""""""""""""""""""""""""""+
%    |                              |
%  3 |                              | 4
%    |                              |
%    |                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   |   
%    +""""""""""""""""""""""""""""""+
%                    2
%

%%%% G E O M E T R Y %%%%
%============================

 %Geometry
xm = 0;  xM = 5;
ym = -1;  yM = 1;

L = xM-xm;  H = yM-ym;
h0 = hh0(jj);  
hH = H/(NNh(jj));
hL = L/ceil(L/hH) ; [hL hH]

x = (xm+hL/2 : hL : xM-hL/2)';  NL = length(x)-1;
y = (ym+hH/2 : hH : yM-hH/2)';  NH = length(y)-1;

% blob size for stokeslet and source dipole
d =dd(jj); dsd = d*1;

% set geometry (position vectors along channel sides)
ex = (xm+hL : hL : xM-hL)';
ey = (ym : hH : yM)';

% we will have points defining all 4 sides
% with forces on all 4 sides
% and sucklets will be placed only on sides 1,2
%
x1 = ex;                y1 =  H/2*ones(size(ex)); 
x2 = ex;                y2 = -H/2*ones(size(ex));
x3 = 0*ones(size(ey));  y3 = ey; 
x4 =  L*ones(size(ey));  y4 = ey;

% where the forces are
xf = [ex;ex;x3;x4];  yf = [y1;y2;y3;y4];

%normals 
nx1 = 0*ex;            ny1 =  ones(size(ex));
nx2 = 0*ex;            ny2 = -ones(size(ex));
nx3 = -ones(size(ey)); ny3 = 0*ey;
nx4 =  ones(size(ey)); ny4 = 0*ey;

n1 = [nx1;nx2;nx3;nx4]; 
n2 = [ny1;ny2;ny3;ny4];

%where the source doublets are
idx1 = find(xf > 0 & xf < L);
idx=[idx1];
xb = [x1;x2]; yb = [y1;y2];

beta0 = 0*xb;   gg=@(z) 1.3*z.^(1) ; 
alpha = hh0(jj)^2 *2.8*NNh(jj)/80*0.9* gg(dsd/dd(1)), bb=0;
% beta0(idx) = -3.4*[ alpha*( ones(size(x1)) ); alpha*( ones(size(x1)))]; 
LL = L*0.5;

beta0(idx) = 3.4*0.95*[ alpha*( ones(size(xb)) ).*((-2*xb).*1/LL^2)]; % linear beta
%beta0(idx) = [ alpha*( ones(size(xb)) ).*((-xb).*3/LL^2)]*2.5; % cubic beta
beta0 = beta0*1;
figure(4),plot(xb,beta0*200)

% KK: why do we redefine xf and yf?
xf = [x1;x2;x3;x4];  yf = [y1;y2;y3;y4];
%                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            
% velocity boundary conditions
%
upf=0*xf;  vpf=0*yf; 
[uex, vex, pex] = permeablechannelexact_semiinf(x4, y4, Da, xM);
upf=[0*x1; 0*x2; 0*x3; uex];
vpf=[0*y1; 0*y2; 0*y3; vex];

U = [upf vpf]; Ulong = U'; Ulong = Ulong(:);
 % compute the forces
xb  = xf(idx);     yb  = yf(idx);     %for sucklets

xe = [x1;x2;x3;x4]; ye = [y1;y2;y3;y4];
wt = [hL*ones(size(x1));hL*ones(size(x2));hH*ones(size(x3));hH*ones(size(x4))];

tmp=[idx*2-1 idx*2]';tmp=tmp(:);  
[Astok,Asuck] = StokSuckCombinedMatrix([xf yf],1,[xf yf],[n1 n2],d,dsd,wt(1:length(xf)),idx,beta0);
Acomb = Astok+Asuck;
% solve for the forces
Ffull = Acomb\Ulong;
f1 = Ffull(1:2:end);  f2 = Ffull(2:2:end);

%fill in the missing bc
[Astok,Asuck] = StokSuckCombinedMatrix([xf yf],0,[xf yf],[[nx1;nx2;nx3] [ny1;ny2;ny3]],d,dsd,wt(1:length(xf)),idx,beta0);
Ufillin = Asuck*Ffull;
uef = [upf;upf(6*NNh:end)];
vef = [vpf;vpf(6*NNh:end)];
uf      = upf + Ufillin(1:2:end);
vf      = vpf + Ufillin(2:2:end);

U = [uf vf]; Ulong = U'; Ulong = Ulong(:);
Ffull = Astok\(Ulong);


% 



% to see the flow due to only the Stokeslet
Ustok = Astok*Ffull;
ustok = Ustok(1:2:end);  vstok = Ustok(2:2:end);  clear Ustok


% ufull = ustok+usuck;  vfull = vstok+vsuck;
% outflows
AA = StokesletVelocityMatrix([xe ye],[xf yf],d,wt(1:length(xf)));
Ufull = AA*Ffull ;  


% Ufake = Astok*Fnew;
ufull = Ufull(1:2:end);  vfull = Ufull(2:2:end);  clear AA

Q1 = sum( ufull(1:NL).*nx1          +vfull(1:NL).*ny1 )*hL;
Q2 = sum( ufull(1+NL:2*NL).*nx2   +vfull(1+NL:2*NL).*ny2)*hL;
Q3 = sum( ufull(1+NL*2:2*NL+NH+2).*nx3+vfull(1+NL*2:2*NL+NH+2).*ny3)*hH;
Q4 = sum( ufull(2*NL+3+NH:end).*nx4 +vfull(2*NL+3+NH:end).*ny4)*hH;

hg = hh0(1)*2;
[xgg,ygg] = meshgrid(0:1*hg:L, -H/2:1*hg:H/2); xg=xgg(:); yg=ygg(:);
ug=0*xgg; vg=ug;
AA = StokesletVelocityMatrix([xg yg],[xf yf],d,wt(1:length(xf)));
Ugrid = AA*Ffull ;  
ug(:) = Ugrid(1:2:end);  vg(:) = Ugrid(2:2:end);  clear AA
speed=sqrt(ug.^2+vg.^2);
if jj==1, u1g = ug; v1g = vg; end
if jj==2, u2g = ug; v2g = vg; end
if jj==3, u3g = ug; v3g = vg; end
if jj==4, u4g = ug; v4g = vg; end
%%
sk=1;
figure(2),%subplot(211)
% pcolor(xgg,ygg,speed),shading interp
plot(xf,yf,'k.'),hold on, hold on,
% quiver(xgg,ygg,ug,vg,0.8,'r','LineWidth',1)
% quiver(xe(1:1*sk:end),ye(1:1*sk:end),ufull(1:1*sk:end)/4,vfull(1:1*sk:end)/4,0,'b','LineWidth',2)
% quiver(xe(1:8:end),ye(1:8:end),usuck(1:8:end),vsuck(1:8:end),0,'r')
% surf(xgg,ygg,-speed),view(2),
hold on
%%%%%%
%%%%%%%%%%%

%%%%%%%%%%%
[upgg,vpgg,ppgg] = permeablechannelexact_semiinf(xgg,ygg, Da, xM);
 
%%%%%%%%%%%%%%%%%%%%%%%%%

hh1=streamline(xgg,ygg,ug ,vg ,xgg(1:end,1 ),ygg( 1:end,1 ));
hh2=streamline(xgg,ygg,ug ,vg ,xgg(1:end,end ),ygg( 1:end,end ));
hh3=streamline(xgg,ygg,-ug ,-vg ,xgg(1:end,end ),ygg( 1:end,end ));

set(hh1,'Color','black');hold off,axis equal,
set(hh2,'Color','black');hold off,axis equal,axis([0-H/2, L, -H/2, H/2  ])
set(hh3,'Color','black');hold off,axis equal,axis([0-H/2, L+H/2, -H*1.5/2, H*1.5/2  ])
     colorbar
hold off  
set(gca, 'FontSize', 16)
pause(0.05)

   figure(2)
hold on

hh1=streamline(xgg,ygg,upgg ,vpgg ,xgg(1:end,1 ),ygg( 1:end,1 ));
hh2=streamline(xgg,ygg,upgg ,vpgg ,xgg(1:end,end ),ygg( 1:end,end ));
hh3=streamline(xgg,ygg,-upgg ,-vpgg ,xgg(1:end,end ),ygg( 1:end,end ));

set(hh1,'Color','red','LineStyle','--');hold off,axis equal,
set(hh2,'Color','blue','LineStyle','--');hold off,axis equal,axis([0-H/2 L -H/2 H/2  ])
set(hh3,'Color','red','LineStyle','--');hold off,axis equal,axis([0-H/2 L+H/2 -H*1.5/2 H*1.5/2  ])
    colorbar
hold off , title('red = exact; black = computed')
set(gca, 'FontSize', 16),pause(0.05)
%%
figure(1)
contourf(xgg,ygg,ppgg)
plot(xf,yf,'k.'),hold on
quiver(xgg, ygg, ug-upgg, vg-vpgg, 1, 'b-')
quiver(xgg, ygg, ug/40, vg/40, 1, 'r-')
hold off,axis equal
% title('sum of both')
format short e
[[Q3 Q1  Q4] sum([Q1 Q2 Q3 Q4])]
netflx(jj)=sum([Q1 Q2 Q3 Q4]);
%end %KK: loop for jj
%%
return
figure(1),clf
loglog(dd,abs(netflx),'-o',dd,abs(netflx(end))*(dd/min(dd)).^1,'--',dd,abs(netflx(end))*(dd/min(dd)),'--','LineWidth',2),
grid,hold on,xlabel('$\epsilon = 2.1 (\Delta s),\ \ \frac{1}{320} < \Delta s < \frac{1}{80}$','interpreter','latex')
ylabel('Net flux'),set(gca,'FontSize',20)
title(['$\beta = -0.0081661 \epsilon$ '],'interpreter','latex')
text(0.015,2.5e-6,'$\epsilon$','interpreter','latex','FontSize',32)
text(0.01,5e-6,'$\epsilon^2$','interpreter','latex','FontSize',32)
%%
return
%%
dd= 0.05/sqrt(5)*[1/3 1/2 1/1.5 1   2 3];
betas=0.0001826*[1/3.1 1/2.013 1/1.504 1 2.039 3.125];
flx= [-9.2505e-04 3.5380e-05 -6.8185e-04 -2.4855e-05  -5.3369e-04 8.6803e-04  ];
plot(dd,betas,'o-',dd, 8.3333e-03*(dd))
%%
% find the normal component of force
fn = f1.*[nx1;nx2;nx3]+f2.*[ny1;ny2;ny3];
%==============================

% now set up a grid and compute the pressure there
m = fix(NL/5);  mm = fix(NH/2);
[xg,yg] = meshgrid(0+xf(m) : hL : L-xf(m), 0-2*mm*hH : 1*hH : H+2*mm*hH);
[ig,jg] = size(xg);

% find the jump for each point of the grid.  For example, if boundary
% point (xf(65),yf(65)) is the closest boundary point to a grid point
% (xg(17,2),yg(17,2)), then jumpidx(17,2) = 65
%
jumpidx = zeros(size(xg));
for k1 = 1 : ig
    for k2 = 1: jg
        tmp = sqrt((xg(k1,k2)-xf(1:2*NL)).^2+(yg(k1,k2)-yf(1:2*NL)).^2);
         jk = find( tmp==min(tmp) , 1 );
         jumpidx(k1,k2) = jk;
    end
end
    
% compute the pressure
pg = stokesletpressure2D(xg,yg,[xf yf],[f1 f2],jumpidx,fn,[n1 n2],H,hL,d,flag);

% change the sign of the one with jumps (not sure why)
if (flag==1), pg = -pg; end
% subtract a constant
%pg = pg - min(min(pg));
%===============================

dpdx = (max(pg(ceil(ig/2),:))-min(pg(ceil(ig/2),:))) ...
      /(max(xg(ceil(ig/2),:)) - min(xg(ceil(ig/2),:)))

% plot
figure(1),clf
mesh(xg,yg,pg),
axis equal,
xlabel('x','FontSize',14),ylabel('y','FontSize',14),zlabel('p','FontSize',14)


figure(2)
plot(xg(ceil(ig/2),:),pg(ceil(ig/2),:),clr)
title('pressure along the channel at y = H/2 ','FontSize',14)
xlabel('x','FontSize',14),ylabel('p','FontSize',14),grid on

figure(3)
plot(yg(:,1),pg(:,1),clr)
title('pressure across the channel at x = 1','FontSize',14)
xlabel('y','FontSize',14),ylabel('p','FontSize',14),grid on



%end%main function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function A = StokesletVelocityMatrix(Xeval,Xforce,del,wt)

Nrow = length(Xeval(:,1));
Ncol = length(Xforce(:,1));
A = zeros(2*Nrow,2*Ncol);

del2 = del*del;
d2 = del2;
mu = 1;

xf = Xforce(:,1);
yf = Xforce(:,2);

x = Xeval(:,1);
y = Xeval(:,2);

for k=1:1:Nrow
   
   rowid = 2*(k-1);    
   rowid1 = rowid+1;   
   rowid2 = rowid+2;
      
   dx = x(k)-xf;
   dy = y(k)-yf;
   r2 = dx.^2+dy.^2;
   R2 = r2+d2;
   
   % H1 = (2*del2)./(r2+del2) - log(r2+del2);
   % H2 = 2./(r2+del2);
   % 
   % blob = 6*d2^2*(2*d2^2-5*d2*r2+r2.^2)./(pi*(r2+d2).^5)
   H1 = (d2*(8*d2^2-4*d2*R2+R2.^2))./(R2).^3 - log(R2);
   H2 = (4*d2^2+d2*R2+R2.^2)./(R2).^3;
   
   %superblob
   % H1 = (2/3*del2*(7*del2^2+r2.^2))./(r2+del2).^3 - log(r2+del2);
   % H2 = 2/3*(15*del2^2+10*del2*r2+3*r2.^2)./(r2+del2).^3;
   
   A(rowid1,1:2:end) = (H2.*dx.*dx + H1 ).*wt;
   A(rowid1,2:2:end) = (H2.*dx.*dy).*wt;
   
   A(rowid2,1:2:end) = (H2.*dy.*dx).*wt;
   A(rowid2,2:2:end) = (H2.*dy.*dy + H1).*wt;
   
end

A = A/(8*pi*mu);
end %StokesletVelocityMatrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%
function pg = stokesletpressure2D(xg,yg,x,f,jumpidx,fn,n,H,hL,d,flag)

%
% (xg,yg) = Mx2 array with points where the pressure will be computed
% x   = Nx2 array where the forces are given
% f   = Nx2 array of forces
% d   = delta 
%
pg = zeros(size(xg));
d2 = d*d;
mu = 1;

% we have to loop over the grid points
[ig,jg] = size(xg);
for k1 = 1:ig
    for k2 = 1:jg
        f1new = f(:,1) - flag*fn(jumpidx(k1,k2))*n(:,1);
        f2new = f(:,2) - flag*fn(jumpidx(k1,k2))*n(:,2);

        dx = xg(k1,k2)-x(:,1);
        dy = yg(k1,k2)-x(:,2);
        r2 = dx.^2 + dy.^2;

        Hs = (r2+2*d2)./(2*pi*(r2+d2).^2);      %Marian's blob
        fdotx = f1new.*dx + f2new.*dy ;
        pg(k1,k2) = sum( fdotx.*Hs/mu )...
               + flag*fn(jumpidx(k1,k2))/hL.*(yg(k1,k2)>0).*(yg(k1,k2)<H);
    end
end


end%stokesletpressure2D
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Ast,Asu] = StokSuckCombinedMatrix(Xeval,choice,Xforce,normals,del,dsd,wt,idx,beta0)

% idx = indices of Xforce where beta is nonzero

Nrow = length(Xeval(:,1));
Ncol = length(Xforce(:,1));
Ast = zeros(2*Nrow,2*Ncol);
Asu = zeros(2*Nrow,2*Ncol);
rhs = zeros(2*Nrow,1);

d2 = del*del;      dsd2 = dsd^2;
mu = 1;

xf = Xforce(:,1);
yf = Xforce(:,2);

n1 = normals(idx,1);
n2 = normals(idx,2);

x = Xeval(:,1);
y = Xeval(:,2);


for k=1:1:Nrow
   
   rowid = 2*(k-1);    
   rowid1 = rowid+1;   
   rowid2 = rowid+2;
      
   dx = x(k)-xf;
   dy = y(k)-yf;
   r2 = dx.^2+dy.^2;
   R2 = r2+d2;
   
   % H1 = (2*d2)./(R2) - log(R2);
   % H2 = 2./(R2);
   S1 = (r2+2*dsd2)./(2*pi*(r2+dsd2).^2); %G'/r
   S2 = -(r2+3*dsd2)./(pi*(r2+dsd2).^3);  %(rG''-G')/r^3

   
   % H1 = (2/3*d2*(7*d2^2+r2.^2))./(8*(r2+d2).^3) - log(r2+d2)/8;
   % H2 = 2/3*(15*d2^2+10*d2*r2+3*r2.^2)./(8*(r2+d2).^3);

   % blob = 6*d2^2*(2*d2^2-5*d2*r2+r2.^2)./(pi*(r2+d2).^5)
   H1 = (d2*(8*d2^2-4*d2*R2+R2.^2))./(R2).^3 - log(R2);
   H2 = (4*d2^2+d2*R2+R2.^2)./(R2).^3;
   
   Ast(rowid1,1:2:end) = (H2.*dx.*dx + H1 ).*wt;
   Ast(rowid1,2:2:end) = (H2.*dx.*dy).*wt;
   
   Ast(rowid2,1:2:end) = (H2.*dy.*dx).*wt;
   Ast(rowid2,2:2:end) = (H2.*dy.*dy + H1).*wt;
   
   %super blob
   m = find(idx==k);
   if isempty(m)==choice
   % if isempty(m)==0
       
   dx = x(k)-xf(idx);
   dy = y(k)-yf(idx);
   r2 = dx.^2+dy.^2;
      ndotx = n1.*dx + n2.*dy;

   %superblob
   % S1 =  (10*dsd2^3 + 5*dsd2^2*r2 + 4*dsd2*r2.^2 + r2.^3)./(2*(r2+dsd2).^4);
   % S2 = -(35*dsd2^3 + 7*dsd2^2*r2 + 5*dsd2*r2.^2 + r2.^3)./(  (r2+dsd2).^5);
   
   Asu(rowid1,idx*2-1) = ((S1(idx).*n1 + S2(idx).*ndotx.*dx).*n1.*beta0(idx).*wt(idx))';
   Asu(rowid1,idx*2  ) = ((0      + S2(idx).*ndotx.*dx).*n2.*beta0(idx).*wt(idx))';
   
   Asu(rowid2,idx*2-1) = ((0      + S2(idx).*ndotx.*dy).*n1.*beta0(idx).*wt(idx))';
   Asu(rowid2,idx*2  ) = ((S1(idx).*n2 + S2(idx).*ndotx.*dy).*n2.*beta0(idx).*wt(idx))';
   end
end

Ast = Ast/(8*pi*mu);
Asu = Asu/(1*mu);
end %function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%
function [usol, vsol, psol] = permeablechannelexact_semiinf(x, y, Da, L_r)
% L_r = L/r in paper

% solve eq (30) in FB
f=@(lam) tan(lam)^2 - tan(lam)/lam + 1 - 2*Da;
lam1 = fzero(f, sqrt(2*Da));

% coefficient from eq (66) of FB
C1 = (-2*L_r)*lam1/sin(lam1);
disp(['lambda_1 = ', num2str(lam1)])

% set up g_1(y) and derivative, eq (40) in FB
g = cot(lam1)*(Da - 0.5)*sin(lam1*y) + 0.5*y.*cos(lam1*y);
g_prime = cot(lam1)*(Da - 0.5)*lam1*cos(lam1*y) + 0.5*cos(lam1*y) - 0.5*lam1*y.*sin(lam1*y);

% ground-state approximate solution from eq (64)-(68), (74)-(75) in FB
psol = C1*exp(lam1*(x - L_r)).*cos(lam1*y);
usol = C1*(-exp(lam1*(x - L_r))).*g_prime/lam1;
vsol = C1*exp(lam1*(x - L_r)).*g;

end