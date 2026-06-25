function [p] = RegStokeslets2D_forcetopressure(y,f,x,ep,blob_num)

% Computes pressures at set of points when given a set of points and
% forces at those points using the Method of Regularized Stokeslets 
% Based on Cortez, SIAM J. Sci Comput. 2001

% Force to velocity code developed by Shilpa Khatri and Ricardo Cortez 
% July 2024 
% Modified by Kristin Kurianski July 26, 2024

%y = (y1,y2) source points of force
%f = (f1,f2) forces at those source points (not force density)
%x = (x1,x2) target points 
%p = pressure evaluated at those target points 
%mu is the viscosity 
%ep is the width of the regularization 

N = size(y,1); %number of source points 
M = size(x,1); %number of target points 

%unpacking the inputs 
y1 = y(:,1);
y2 = y(:,2);
f1 = f(:,1);
f2 = f(:,2); 
x1 = x(:,1); 
x2 = x(:,2); 

%initializing the velocity 
p = zeros(M,1);

%loop over source points    
for k = 1:N 

    %distance between target and source points 
    XY1 = x1(:) - y1(k); 
    XY2 = x2(:) - y2(k);  
    R2 = XY1.^2 + XY2.^2 + ep^2; 
    R = sqrt( R2 ); 

    %computing the pressure
    %Note: S1 = G'(r)/r
    %p(r) = (f . x)S1
    [~, ~, S1, ~, ~] = reg_fncs_withdoublet(ep,R,blob_num);

    fdotXY = f1(k)*XY1 + f2(k)*XY2; 
    
    p(:) = p(:) + fdotXY.*S1;

end

end