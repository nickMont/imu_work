function [xk,RR2,L2] = updateRandL(xk,RR,L1)
RR2 = euler2dcm(xk(7:9))*RR;
xk(7:9)=zeros(3,1);
L2 = L1+xk(16:18);
xk(16:18)=zeros(3,1);
end

