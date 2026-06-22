function [chi_inside] = IBDLcharacteristic(xmin, ymin, Lx, Ly,Nx, Ny,...
    X0,  unitnormals , ds_s )

dx=Lx/Nx; 
dy=Ly/Ny;

xg=dx*(0:Nx-1)+xmin;
yg=dy*(0:Ny-1)+ymin;
[xg,yg]=ndgrid(xg,yg);



S1=spreadmatrix_vc_vec(X0, dx, Nx,Ny,xmin,ymin);
Spread1=S1.*ds_s'/dx^2;
Sn1=Spread1*unitnormals;
[Dx, Dy]=dx2d(Nx,Ny);
Dx=1/(2*dx)*Dx;
Dy=1/(2*dx)*Dy;
DivSn1=Dx*Sn1(:,1)+Dy*Sn1(:,2);
DivSn1=reshape(DivSn1, Nx, Ny,1);
chi1=helmholtz_solve_FD(DivSn1,1,0, Lx,Ly, dx,dy);
chi1=chi1-chi1(1,1);
interiorofdiscreteshape1=find(chi1(:)>=1/2);
chi_inside=zeros(size(xg));
chi_inside(interiorofdiscreteshape1)=1;



end