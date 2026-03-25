function [T02,T04,T06,T12,T14,T16,T24,T26] = sd_segment_terms(x,y,ep)
%
%   * x = [x1, x2] are the evaluation/target points
%   * y = [yj,yk] are the segment end points
%   * ep is the epsilon value for the segment
% 
% Inputs x and ycan be ( )x2 matrices, so that you can output T's for
% multiple segments and evaluation points. ep should be the same length as
% y.
%
% Developed by Michaela Kubacki and Ricardo Cortez March 2026

yj = y(1,:);
yk = y(2,:);
x1 = x(:,1);
x2 = x(:,2);

xj1 = x1 - yj(1);
xj2 = x2 - yj(2);
xk1 = x1 - yk(1);
xk2 = x2 - yk(2);

x0 = [xj1;xj2];
x1 = [xk1;xk2];

ep2 = ep.^2;

R0 = sqrt(xj1.^2 + xj2.^2 + ep2);
R1 = sqrt(xk1.^2 + xk2.^2 + ep2);

v = yj-yk;
L = norm(v);

xjdotv = xj1*v(1) + xj2*v(2);
xkdotv = xk1*v(1) + xk2*v(2);

Log0 = log(R0);
Log1 = log(R1);

c = R0.^2 - (xjdotv).^2/L^2;
csqrt = sqrt(c);

% xi = aL + (xjdotv)/L
xi1 = L + (xjdotv)/L; % when a = 1
xi0 = xjdotv/L; % when a = 0

Atan1 = atan2(xi1,csqrt);
Atan0 = atan2(xi0,csqrt);

R02 = R0.^2; R04 = R02.*R02; 
R12 = R1.^2; R14 = R12.*R12; 

% T02 calculated directly
T02 = (Atan1-Atan0)./(L*csqrt);

% Rest of T0's calculated using recursion
T04 = (0.5*T02 + 0.5*xi1./(R1.^2)/L - 0.5*xi0./(R0.^2)/L)./c;

T06 = (0.75*T04 + 0.25*xi1./(R1.^2).^2/L - 0.25*xi0./(R0.^2).^2/L )./c;

% Rest of T's Calculated Using the recursion
T12 = (Log1 - Log0)/L^2 - xjdotv.*T02/L^2;

T14 = -0.5/L^2 * (1./R12 - 1./R02) - xjdotv.*T04/L^2;
T24 = -0.5/L^2 * (1./R12) + 0.5/L^2 * T02 - xjdotv.*T14/L^2;

T16 = -0.25/L^2 * (1./R14 - 1./R04) - xjdotv.*T06/L^2;
T26 = -0.25/L^2 * (1./R14) + 0.25/L^2 * T04 - xjdotv.*T16/L^2;

