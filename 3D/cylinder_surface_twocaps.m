function [Px,Py,Pz,dArea,Ntube,Ncap]=cylinder_surface_twocaps(Nt,L,a)

%discretization of circle
dt = 2*pi/Nt;
t = 0:dt:2*pi-dt/2;
t = t';

%using code from RC to discretize cap of circle 
%results in discretization with dArealid approx (a*dt)^2*pi/3 
%h is the actual discretization computed from optimal a*dt and dArealid = h^2*pi/3
[xlid,ylid,zlid,dArealid,h] = make_lid(a,a*dt); 
Ncap = length(xlid);
Nloop = sqrt(Ncap/3); %number of annuli in discretization of lid 
Nouter = (2*Nloop-1)*3; %number of points on outer most annulus 

%discretization of length of cylinder - note dL is computed to be approx h 
dL = h;
NL = ceil(L/dL); 
dL = L/NL;
s = dL/2:dL:L; 
s = s'; 


counter = 1; 
%points on the tube of cylinder 
for i=1:length(s)
    for j = 1:Nouter %set the number of points on the radial direction to be the same as on lid 
        y1(counter,1) = s(i); 
        y2(counter,1) = ylid(end-Nouter+j);
        y3(counter,1) = zlid(end-Nouter+j);
        counter = counter + 1; 
    end
end

%area of cell on tube 
da = a*2*pi*dL/Nouter;

Ntube = length(y1); 
dArea(1:Ntube,1) = da; 



%% using builtin pde toolbox to discretize cap of circle 
% xcir = a*cos(t); 
% ycir = a*sin(t); 
% pcircle = polyshape(xcir,ycir);
% model = createpde; 
% tr = triangulation(pcircle);
% cnodes = tr.Points'; 
% celements = tr.ConnectivityList'; 
% geom = geometryFromMesh(model,cnodes,celements); 
% FEMesh = generateMesh(model,"Hmin",dt,GeometricOrder="linear");
% meshNodes = FEMesh.Nodes';
% meshTris = FEMesh.Elements';
% 
% %finding center of triangles in cap mesh 
% for i = 1:length(meshTris)
%     
%     ind1 = meshTris(i,1); 
%     ind2 = meshTris(i,2); 
%     ind3 = meshTris(i,3);
% 
%     tri_center(i,:) = meshNodes(ind1,:) + meshNodes(ind2,:) + meshNodes(ind3,:);
%     tri_center(i,:) = tri_center(i,:)/3;
% 
% end
% Ncap = length(tri_center); 
% 
% %area of triangles 
% T = meshTris; 
% tri_dArea = 0.5*( (meshNodes(T(:,2),1)-meshNodes(T(:,1),1)).*(meshNodes(T(:,3),2)-meshNodes(T(:,1),2)) - ...
%     (meshNodes(T(:,3),1)-meshNodes(T(:,1),1)).*(meshNodes(T(:,2),2)-meshNodes(T(:,1),2)) );
% %----------------------------------------------------

%putting together tube + 2 caps 
Px = [y1;0*ones(Ncap,1);L*ones(Ncap,1)]; 
Py = [y2;ylid;ylid]; 
Pz = [y3;zlid;zlid]; 
dArea = [dArea;dArealid;dArealid];





