function [xkp1,Pkp1,RBI_out] = ekfUnknownLever(x0,Pk,imuMeas,gpsMeas,systemParams,RBI0)
%get all time stamps. GPS is assumed to be the last measurement in terms of
%time
%Sk and nuj are only used in MMKF

%gpsMeas is [t0;gps1;gps2];

maxLeverPossible=[1;1;1];  %change to 1e7 to ignore

Sk=[]; nuj=[]; xkp1=x0; Pkp1=Pk; RBI_out=RBI0;

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
Qimu=systemParams.Rimu;

L_s2p=systemParams.Ls2p;
L_cg2p=systemParams.Lcg2p;
tauA=systemParams.tA;
tauG = systemParams.tG;

%Limu=[0;0;0];

if numImuMeas>1
    eye3=eye(3); zer3=zeros(3,3);
    
    Gt=zeros(18,12);
    Ft=eye(18);
    xk = x0;
    
    %Run CF
    RR=RBI0;
    dtsum=0;
    [xk,RR]=updateRBI(xk,RR);
    for i=1:numImuMeas
        dt=timeVec(i+1)-timeVec(i); %go through all imu and then gpstoimu time
        dtsum=dtsum+dt;
        %[exvCF,Pin]=complimentaryFilterWithLever(dt,exvCF,imuMeas{i},Pin,biasEstimate,L_ab0-biasState(16:18));
        %[exvCF,Pin]=complimentaryFilterWithLever_v2(dt,RR,exvCF,imuMeas{i},Pin,biasEstimate,L_ab0+biasState(16:18));
        fB = imuMeas{i}(5:7);
        wB = imuMeas{i}(2:4);
        
        xk = f_imu_dyn_unknownLever(dt,xk,RR,fB,wB,zeros(12,1),tauA,tauG);
        F_local = complexStep(@(xVar) f_imu_dyn_unknownLever(dt,xVar,RR,fB,wB,zeros(12,1),tauA,tauG),xk,1e-10);
        Ft=F_local*Ft;
        Gt = Gt + complexStep(@(vVar) f_imu_dyn_unknownLever(dt,xk,RR,fB,wB,vVar,tauA,tauG),zeros(12,1),1e-10);
        [xk,RR]=updateRBI(xk,RR);
        %Limu=[0;0;0];
    end
    
    xbar = xk;
    zk=[gpsMeas{1}(2:4); unit3(gpsMeas{1}(5:7))];  %pose meas in local frame
    Hk = complexStep(@(xVar) h_imu_meas(xVar,RR,zeros(6,1),L_cg2p,L_s2p),xbar,1e-10);
    Rk=.0002*[eye3 zer3; zer3 eye3];
    Pbar=Ft*Pk*Ft'+Gt*Qimu*Gt';
    %diagFPF = diag(Ft*Pk*Ft');
    %diagGQG = diag(Gt*Qimu*Gt');
    %diagP=diag(Pbar);
    Sk=Rk+Hk*Pbar*Hk';
    Wk=Pbar*Hk'*inv(Sk);
    nuj=zk-h_imu_meas(xbar,RR,zeros(6,1),L_cg2p,L_s2p) 
    xkp1=xbar+Wk*nuj;
    
    Pkp1 = (eye(18)-Wk*Hk)*Pbar*(eye(18)-Wk*Hk)' + Wk*Rk*Wk'; 
    [xkp1,RBI_out]=updateRBI(xk,RR);
end

end

