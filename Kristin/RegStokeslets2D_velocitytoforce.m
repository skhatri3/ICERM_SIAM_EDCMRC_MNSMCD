function [f] = RegStokeslets2D_velocitytoforce(y,x,u,ep,mu,blob_num)

% Computes forces at set of points when given a set of points and
% velocities at those points using the Method of Regularized Stokeslets 
% Based on Cortez, SIAM J. Sci Comput. 2001

% Developed by Shilpa Khatri and Ricardo Cortez 
% July 2024 

%y = (y1,y2) source points
%f = (f1,f2) forces at those source points (not force density) 
%x = (x1,x2) target points 
%u = (u1,u2) velocity at those target points 
%mu is the viscosity 
%ep is the width of the regularization 

%ACCURACTY NOTE: If inf norm of condition number greater than 10^r (r>0)
%                then losing r digits of accuracy in matrix solve. Warning
%                message will display if r >= 8.

N = size(y,1); %number of source points 
M = size(x,1); %number of target points 

%unpacking the inputs 
y1 = y(:,1);
y2 = y(:,2);
u1 = u(:,1);
u2 = u(:,2); 
x1 = x(:,1); 
x2 = x(:,2); 

%initializing the forces
f1 = zeros(N,1);
f2 = zeros(N,1); 

%building matrix 


%loop over source points    
for k = 1:N 

    %distance between target and source points 
    XY1(:,k) = x1(:) - y1(k); 
    XY2(:,k) = x2(:) - y2(k);  

end

R2 = XY1.^2 + XY2.^2 + ep^2; 
R = sqrt( R2 ); 

[H1, H2] = reg_fncs(ep, R, blob_num); 

M11 = H1 + H2.*XY1.*XY1; 
M22 = H1 + H2.*XY2.*XY2; 
M12 = H2.*XY1.*XY2; 
%M21 = H2.*XY2.*XY1; 

M = [M11 M12; M12 M22]/4/pi/mu;

%solving for force when given a velocity 
uu = [u1;u2];
disp(['log10 of condition number of M: ', num2str(log10(condest(M)))])
disp(['condition number of M: ', num2str(condest(M))])
disp(['size of M: ', num2str(size(M))])
log10cond = log10(condest(M));
if log10cond >= 8
    disp(['Warning: High condition number. Losing roughly ', num2str(log10cond), ' digits of accuracy.'])
end
ff = M\uu; 
f1 = ff(1:N); 
f2 = ff(N+1:end); 

%repacking output 
f = [f1,f2]; 

