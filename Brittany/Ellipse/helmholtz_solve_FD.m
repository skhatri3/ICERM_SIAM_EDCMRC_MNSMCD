%
%  Solve the helmholtz equation
%    a(u_{xx} + u_{xx}) - b u + f = 0
%  in a periodic box of size (Lx,Ly)
% using eigenvalues for a FD matrix instead 
%
function u=helmholtz_solve_FD(f,a,b,Lx,Ly, dx, dy)
    
    % record the number of grid points in each direction
    %
    szf=size(f);
    nx=szf(1);  
    ny=szf(2);
      
    % compute the wave numbers
    %
    N1x =  floor((nx-1)/2);
    N2x = (nx/2)*ones(rem(nx+1,2));
    freqx =(2*pi/Lx)* [(0:N1x)  N2x (-N1x:-1)]';

    N1y =  floor((ny-1)/2);
    N2y = (ny/2)*ones(rem(ny+1,2));
    freqy = (2*pi/Ly)*[(0:N1y)  N2y (-N1y:-1)]';

    [k1 k2]=ndgrid(freqx,freqy);

    % compute eigenvalues of the operator
    %
    kk = a/dx^2*(2*cos(k1*dx) + 2*cos(k2*dy)-4) - b;
        
    % if b=0, the problem is singular, adjust this eigenvalue
    %
    if( b==0 )
      kk(1) = 1;
    end
    
    % transform f, and project off mean if problem is singular
    %
    fhat = fft2(f);
    if( b==0 )
      fhat(1) = 0;
    end
    
    % do the solve
    %
    uhat = -fhat./kk;  % note the minus sign because I'm solving -L=f
    u    = real( ifft2(uhat) );
    
