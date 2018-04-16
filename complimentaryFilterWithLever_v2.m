function [gammaxv2,Pin] = complimentaryFilterWithLever_v2(dt,R0,gammaxv1,imuMeas,Pin,biasEstimate,L_ab)
%biasEstimate is [ba;bg]
A=eye(6); A(1:3,4:6)=dt*eye(3);
B=[.5*dt^2*eye(3);dt*eye(3)];
gammaxv2=zeros(9,1);
wB=imuMeas(2:4)-biasEstimate(4:6);
Rhat=R0*(eye(3)+hatmat(gammaxv1(1:3)));
gammaxv2(1:3)=gammaxv1(1:3) + dt*eye(3)*wB;
%emid=(exv2(1:3)+exv1(1:3)); %appears to be more accurate to use e(k) than to average e(k),e(k+1)
%R_body2local=calculate_R_from_euler(emid)';
g=[0;0;-9.81];
%g=[0;0;0];
gammaxv2(4:9)=A*gammaxv1(4:9) + B*Rhat'*(imuMeas(5:7)-biasEstimate(1:3)) + ...
    B*g - B*cross(wB,cross(wB,L_ab));  
%exv2(1:3)=exv1(1:3)+dt*eye(3)*imuMeas(2:4);
end

