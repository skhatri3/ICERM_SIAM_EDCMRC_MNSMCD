function [Q0,Q1,T02,T12,T22,T32] = st_segment_terms(x,y,ep)
%
%   * x = [x1, x2] are the evaluation/target points
%   * y = [yj,yk] are the segment end points
%   * ep is the epsilon value for the segment
% 
% Inputs x and ycan be ( )x2 matrices, so that you can output T's Q's for
% multiple segments and evaluation points. ep should be the same length as
% y.
%
% Developed by Michaela Kubacki and Ricardo Cortez March 2026

yj = y(1,:); % left endpoint
yk = y(2,:); % right endpoint

x1 = x(:,1);
x2 = x(:,2);

xj1 = x1 - yj(1);
xj2 = x2 - yj(2);
xk1 = x1 - yk(1);
xk2 = x2 - yk(2);

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

% Q0 calculated directly
Q0=1/L^2*( xkdotv.*(Log1-1) - xjdotv.*(Log0-1) + L*csqrt.*(Atan1-Atan0) );

% T02 calculated directly
T02 = (Atan1-Atan0)./(L*csqrt);

% Using hand-derived formula (verified) 
Q1 = (2*R1.^2.*Log1 - R1.^2)/(4*L^2) - (2*R0.^2.*Log0 - R0.^2)/(4*L^2) - 1/(L^2)*xjdotv.*Q0;

% Use recursion formulas for remaining terms
T12 = (Log1-Log0)/L^2 - xjdotv.*T02/L^2;
T22 = Log1/L^2 - Q0/L^2 - xjdotv.*T12/L^2;
T32 = Log1/L^2 - 2*Q1/L^2 - xjdotv.*T22/L^2;