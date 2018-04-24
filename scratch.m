figure(1);clf;
pli=1:1:length(leverEst);
hold on
plot(thist0(1+pli),leverEst(1,pli)-Lab_true(1),':r')
hold on
plot(thist0(1+pli),leverEst(2,pli)-Lab_true(2),'--b')
hold on
plot(thist0(1+pli),leverEst(3,pli)-Lab_true(3),'-g')
hold on
legend('\delta l_x','\delta l_y','\delta l_z')
axis([0 20 -1.5 1.5])
xlabel('Time (s)')
ylabel('Error (m)')
title('Lever arm estimation error with IMMEKF (10cm discretization)')
grid on
figset