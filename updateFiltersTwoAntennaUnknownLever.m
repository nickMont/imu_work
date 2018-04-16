function [ewxvUp,biasState,Pout,Sk,nuj] = updateFiltersTwoAntennaUnknownLever(dtt,ewxv,Pin,imuMeas,gpsMeas,biasState,L_ab0,L21)
%get all time stamps. GPS is assumed to be the last measurement in terms of
%time
%Sk and nuj are only used in MMKF

%gpsMeas is [t0;gps1;gps2];

maxLeverPossible=[1;1;1];  %change to 1e7 to ignore

Sk=[]; nuj=[];

%xState = [dP;dV;eul;bg;ba;dL];
numImuMeas=length(imuMeas);
timeVec=zeros(numImuMeas+1,1);
if numImuMeas>=1
    for i=1:numImuMeas
        timeVec(i)=imuMeas{i}(1);
    end
end
timeVec(end)=gpsMeas{1}(1);

Re=637800;

biasEstimate=biasState(10:15);
%biasEstimate=zeros(6,1);
Qimu=diag([10^-7*ones(3,1); 10^-8*ones(3,1)]);

if numImuMeas>1
    RR=calculate_R_from_euler(ewxv(1:3));
    %exvCF=[ewxv(1:3);ewxv(7:12)];
    exvCF=[zeros(3,1);ewxv(7:12)];
    eye3=eye(3); zer3=zeros(3,3);
    
    %Run CF
    for i=1:numImuMeas
        dt=timeVec(i+1)-timeVec(i); %go through all imu and then gpstoimu time
        %[exvCF,Pin]=complimentaryFilterWithLever(dt,exvCF,imuMeas{i},Pin,biasEstimate,L_ab0-biasState(16:18));
        [exvCF,Pin]=complimentaryFilterWithLever_v2(dt,RR,exvCF,imuMeas{i},Pin,biasEstimate,L_ab0+biasState(16:18));
        Ak=[zer3 eye3 zer3 zer3 zer3 zer3
            zer3 zer3 -RR*hatmat(imuMeas{i}(5:7)) zer3 RR zer3
            zer3 zer3 -hatmat(imuMeas{i}(2:4)) eye3 zer3 zer3
            zeros(9,18)];
        Ft=(eye(18)+dt*Ak)*Ft;
    end
    
    %note: biasState = [de, dx, dv, ba, bg]
    %Gammak: note that accel noise enters in body, not world, frame
    Gt=[zer3 dtt*eye3
        RR*dtt^2*eye3/2 zer3
        RR*dtt*eye3 zer3
        eye3 zer3
        zer3 eye3
        zer3 zer3];
    Fb = hatmat(imuMeas{numImuMeas}(5:7));
    Omegab=hatmat(imuMeas{numImuMeas}(2:4));
    L1b=hatmat(L_ab0);
    Fk23 = -RR*Fb;
    Fk33 = -Omegab;
    Ak=[zer3 eye3 zer3 zer3 zer3 zer3
        zer3 zer3 Fk23 zer3 RR zer3
        zer3 zer3 Fk33 eye3 zer3 zer3
        zeros(9,18)];
    Fk=eye(18)+dtt*Ak;
    zk=[(exvCF(4:6)-gpsMeas{1}(2:4)); (gpsMeas{1}(2:4)-gpsMeas{1}(5:7))];  %pose meas in local frame
    Rk=.02*[eye3 zer3; zer3 eye3];
    Hk=[eye3 zer3 -RR*L1b zer3 zer3 RR
        zer3 zer3 -RR*hatmat(L21) zer3 zer3 zer3];
    biaskp1=Fk*biasState;
        %biaskp1(6)=biaskp1(6)+9.81*dtt^2/2; biaskp1(9)=biaskp1(9)+9.81*dtt;
    Pkp1=Fk*Pin*Fk'+Gt*Qimu*Gt';
    Sk=Rk+Hk*Pkp1*Hk';
    Kg=Pkp1*Hk'*inv(Sk);
    nuj=zk-Hk*biaskp1;
    biasState=biaskp1+Kg*nuj;
    
    %TESTING: Saturation limit for known bounds
    biasState(16:18) = vectorSaturationF(biasState(16:18),maxLeverPossible);
    
end

if numImuMeas>1 %split for initialization
    ewxvUp=[exvCF(1:3);imuMeas{numImuMeas}(2:4);exvCF(4:9)];
    %ewxvUp=[exvCF(1:3);imuMeas{numImuMeas}(2:4);gpsMeas{1}(5:7)-biasState(1:3);exvCF(7:9)];
    Pout=Pkp1-Kg*Hk*Pkp1;
    %ewxvUp=[biasState(7:9)+ewxv;imuMeas{numImuMeas}(2:4);biasState(1:6)+defZeroLoc(4:9)];
else
    ewxvUp=ewxv;
    Pout=Pin;
end

end

