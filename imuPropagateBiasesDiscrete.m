function [b1b2e1e2_update, accelBias, gyroBias] = imuPropagateBiasesDiscrete(b1b2e1e2,dt,constValsVec,noiseVals)
%noiseVals is a 9x1 noise vector for bias1,bias2,epsilon2
%constValsVec is a 9x1 time constant vector
%Biases propagated via first-order approximation of ODE

b1b2e1e2_update=zeros(12,1);
b1b2e1e2_update(1:6)=b1b2e1e2(1:6) - inv(diag(constValsVec(1:6)))*b1b2e1e2(1:6)*dt + noiseVals(1:6);
b1b2e1e2_update(7:9)=b1b2e1e2(7:9);
b1b2e1e2_update(10:12)=b1b2e1e2(10:12) - inv(diag(constValsVec(7:9)))*b1b2e1e2(10:12)*dt + noiseVals(7:9);
accelBias=b1b2e1e2_update(1:3) + b1b2e1e2_update(4:6);
gyroBias=b1b2e1e2_update(7:9) + b1b2e1e2_update(10:12);

end

