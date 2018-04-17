function [b1b2e1e2_update, accelBias, gyroBias] = imuPropagateBiasesDiscreteOnebias(b1b2e1e2,dt,constValsVec,noiseVals)
%noiseVals is a 9x1 noise vector for bias1,bias2,epsilon2
%constValsVec is a 9x1 time constant vector
%Biases propagated via first-order approximation of ODE

b1b2e1e2_update=zeros(6,1);
b1b2e1e2_update(1:3)=b1b2e1e2(1:3) - inv(diag(constValsVec(1:3)))*b1b2e1e2(1:3)*dt + noiseVals(1:3);
b1b2e1e2_update(4:6)=b1b2e1e2(4:6) - inv(diag(constValsVec(4:6)))*b1b2e1e2(4:6)*dt + noiseVals(4:6);
accelBias=b1b2e1e2_update(1:3);
gyroBias=b1b2e1e2_update(4:6);

end

