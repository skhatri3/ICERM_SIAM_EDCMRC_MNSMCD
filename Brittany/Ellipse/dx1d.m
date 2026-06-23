function [Dx] = dx1d(nx)
%2D centered difference approx to Dx with periodic boundary conditions
%without 1/(2h) factor

%size is nx ny    x    nxny

ex=ones(nx,1);
Dx=spdiags([-1*ex ex], [-1 1], nx,nx);
Dx(1,end)=-1;
Dx(end,1)=1;





end

