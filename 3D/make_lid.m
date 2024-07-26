function [xl,yl,zl,dArea,h] = make_lid(R,h0)

% Discretize a circular disc of radius R using a
% desired discretization size of h0 (1D)  

Nloop = ceil(R/h0); %number of annuli
h = R/Nloop;        %actual h (rather than desired hO)    

Nlid = 3*Nloop^2;

xl = zeros(Nlid,1);
yl = xl; zl = xl;

yc = 0; zc = 0; %coordinate of center of circle, assume orgin 

ncnt = 0;

for k=0:Nloop-1
   NThisLoop = 3*(2*k+1);
   RThisLoop = (k+1/2)*h;
   dtheta = 2*pi/NThisLoop;
   for j=0:(NThisLoop-1)
     ncnt=ncnt+1;
     yl(ncnt)=RThisLoop*cos(j*dtheta) + yc;
     zl(ncnt)=RThisLoop*sin(j*dtheta) + zc;
     xl(ncnt)=0;
   end
end

dArea = h^2*pi/3*ones(ncnt,1);

%number of points on the outermost loop is 3*(2*Nloop + 1) 