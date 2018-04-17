function [xkp1,Pkp1,RBI_out,Limu] = ukfUnknownLever(x0,Pk,imuMeas,gpsMeas,L0,systemParams,RBI0)
%get all time stamps. GPS is assumed to be the last measurement in terms of
%time
%Sk and nuj are only used in MMKF

%gpsMeas is [t0;gps1;gps2];

maxLeverPossible=[1;1;1];  %change to 1e7 to ignore

Sk=[]; nuj=[]; xkp1=x0; Pkp1=Pk; Limu=L0; RBI_out=RBI0;

%xState = [dP;dV;eul;bg;ba;dL];
numImuMeas=length(imuMeas);
timeVec=zeros(numImuMeas+1,1);
if numImuMeas>=1
    for i=1:numImuMeas
        timeVec(i)=imuMeas{i}(1);
    end
end
timeVec(end)=gpsMeas{1}(1);
%biasEstimate=zeros(6,1);

L_s2p=systemParams.Ls2p;
L_cg2p=systemParams.Lcg2p;
tauA=systemParams.tA;
tauG = systemParams.tG;
Qimu = systemParams.Rimu;
Limu=L0;

if numImuMeas>1
    eye3=eye(3); zer3=zeros(3,3);
    
    xk = x0;
    %Run CF
    RR=RBI0;
    [xk,RR,Limu]=updateRandL(xk,RR,Limu);
    %Limu=[0;0;0];
    for i=1:numImuMeas
        dt=timeVec(i+1)-timeVec(i); %go through all imu and then gpstoimu time
        fB = imuMeas{i}(5:7);
        wB = imuMeas{i}(2:4);
        
        [xk,Pk]=ukfPropagate(dt,xk,Pk,Qimu, RR,fB,wB,L0,tauA,tauG);
        [xk,RR,Limu]=updateRandL(xk,RR,Limu);
        %Limu=[0;0;0];
    end
    
    xbar=xk;
    Pbar=Pk;
    
    %note: biasState = [de, dx, dv, ba, bg]
    %Gammak: note that accel noise enters in body, not world, frame
    
    zk=[gpsMeas{1}(2:4); unit3(gpsMeas{1}(5:7)-gpsMeas{1}(2:4))];  %pose meas in local frame
    Rk=.0002*[eye3 zer3; zer3 eye3];
    [zbar,Pxz,Pzz]=ukfMeasure(xbar,Pbar,Rk,RR,L_cg2p,L_s2p);
    nuj = zk-zbar
    xkp1 = xbar+Pxz*inv(Pzz)*nuj;
    Pkp1 = Pbar - Pxz*inv(Pzz)*Pxz';
    
    %TESTING: Saturation limit for known bounds
    %biasState(16:18) = vectorSaturationF(biasState(16:18),maxLeverPossible);
    [xkp1,RBI_out,Limu]=updateRandL(xkp1,RR,Limu);
end


end

