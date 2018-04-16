function xvdot = scratch_ode(t,xv)
global u;
xvdot=[0;0];
xvdot(1)=xv(2);
xvdot(2)=u-7.35;


end

