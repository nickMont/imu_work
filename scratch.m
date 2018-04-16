clear;clc;

xv=[0;.1];
global u;
k=.3;
kd=.2;
tstore=[]; xvstore=[];
tsamp=1/15;
for ij=1:151
   xvdes=[0;0]; 
   u=k*(xvdes(1)-xv(1))+kd*(xvdes(2)-xv(2))+7.35;
   [tsim,xvsim]=ode45('scratch_ode',tsamp*(ij-1)+[0 tsamp],xv);
   xv=xvsim(end,:)';
   tstore=[tstore;tsim];
   xvstore=[xvstore;xvsim];
end

figure(1);clf;
plot(tstore,xvstore(:,1))
axis([0 10 -.25 .25])

