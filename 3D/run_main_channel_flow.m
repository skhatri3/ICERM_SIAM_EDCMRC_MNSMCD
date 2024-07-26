%run convergence study on main_channel_flow 

clear all
close all

% Plotting figures 
skip = 8; %for quiver plots 
set(0,'defaultaxesfontsize',20,'defaultaxeslinewidth',2.0,...
      'defaultlinelinewidth',2.0,'defaultlinemarkersize',10.0)

%% vary resolution 

Nt = [5,10,20,40,80,160]; 
c_ep = 0.1; 

for i = 1:length(Nt)

    [l2error1(i),l2error2(i),l2error3(i),maxerror1(i),maxerror2(i),maxerror3(i)]=main_channel_flow(Nt(i),c_ep);

end

figure(100) 
loglog(Nt,l2error1,'k','LineWidth',2)
hold on 
loglog(Nt,l2error2,'b')
loglog(Nt,l2error3,'r')
loglog(Nt,maxerror1,'k--','LineWidth',2)
loglog(Nt,maxerror2,'b--')
loglog(Nt,maxerror3,'r--')
title('Convergence in Resolution')
legend('l2 x','l2 y','l2 z','max x','max y','max z')

%% vary blob width 
Nt = 40; 
c_ep = [0.1:0.1:2];

for i = 1:length(c_ep)

    [l2error1ep(i),l2error2ep(i),l2error3ep(i),maxerror1ep(i),maxerror2ep(i),maxerror3ep(i)]=main_channel_flow(Nt,c_ep(i));

end


figure(101) 
loglog(c_ep,l2error1ep,'k','LineWidth',2)
hold on 
loglog(c_ep,l2error2ep,'b')
loglog(c_ep,l2error3ep,'r')
loglog(c_ep,maxerror1ep,'k--','LineWidth',2)
loglog(c_ep,maxerror2ep,'b--')
loglog(c_ep,maxerror3ep,'r--')
title('Convergence in Blob Size')
legend('l2 x','l2 y','l2 z','max x','max y','max z')







