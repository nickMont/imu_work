function [xk,RR2] = updateRBI(xk,RR)
RR2 = euler2dcm(xk(7:9))*RR;
xk(7:9)=zeros(3,1);
end

