function [exv2,Pin] = complimentaryFilterWithLever(dt,exv1,imuMeas,Pin,biasEstimate,L_ab)
%biasEstimate is [ba;bg]
A=eye(6); A(1:3,4:6)=dt*eye(3);
B=[.5*dt^2*eye(3);dt*eye(3)];
R_body2local=calculate_R_from_euler(exv1(1:3))';
%yes, this needs the transpose as rpy are the local2body angles
exv2=zeros(9,1);
wB=imuMeas(2:4)-biasEstimate(4:6);
exv2(1:3)=exv1(1:3) + dt*eye(3)*R_body2local*wB;
%emid=(exv2(1:3)+exv1(1:3)); %appears to be more accurate to use e(k) than to average e(k),e(k+1)
%R_body2local=calculate_R_from_euler(emid)';
g=[0;0;-9.81];
%g=[0;0;0];
exv2(4:9)=A*exv1(4:9) + B*R_body2local*(imuMeas(5:7)-biasEstimate(1:3)) + ...
    -B*g - B*R_body2local*cross(R_body2local'*wB,L_ab);  
%exv2(1:3)=exv1(1:3)+dt*eye(3)*imuMeas(2:4);
end

