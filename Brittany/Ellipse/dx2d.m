function [Dx, Dy] = dx2d(nx,ny)
%2D centered difference approx to Dx with periodic boundary conditions
%without 1/(2h) factor

%size is nx ny    x    nxny

ex=ones(nx,1);
dx=spdiags([-1*ex ex], [-1 1], nx,nx);
dx(1,end)=-1;
dx(end,1)=1;

ey=ones(ny,1);
dy=spdiags([-1*ey ey], [-1 1], ny,ny);
dy(1,end)=-1;
dy(end,1)=1;

Ix=speye(nx);
Iy=speye(ny);

Dx=kron(Iy,dx);
Dy=kron(dy,Ix);



end

