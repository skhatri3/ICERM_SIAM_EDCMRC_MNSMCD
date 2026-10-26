%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ug,vg,pg] = permeablechannelexact_FB(x,y,params)

eta = params.eta;
Da = params.Da;

%%%%%%%%%%%
g  = cot(eta)*(Da-0.5)   *sin(eta*y)+0.5*y.*cos(eta*y);
gp = cot(eta)*(Da-0.5)*eta*cos(eta*y)+0.5*cos(eta*y)-0.5*eta*y.*sin(eta*y);

C = -1;   
pg = C*sinh(eta*x) .* cos(eta*y);
ug = -C/eta * cosh(eta*x) .* gp;
vg =  C    * sinh(eta*x) .* g;
%%%%%%%%%%%%%%%%%%
end