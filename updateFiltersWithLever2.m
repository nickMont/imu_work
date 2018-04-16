function [ewxvUp,biasState,Pout,Sk,nuj] = updateFiltersWithLever2(dtt,ewxv,Pin,imuMeas,gpsMeas,biasState,L_ab)
%ewxv represents state about which to linearize
%get all time stamps. GPS is assumed to be the last measurement in terms of
%time
%updateFiltersWithLever2 does not incorporate attitude measurements
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
    exvCF=[ewxv(1:3);ewxv(7:12)];
    for i=1:numImuMeas
        dt=timeVec(i+1)-timeVec(i); %go through all imu and then gpstoimu time
        [exvCF,Pin]=complimentaryFilterWithLever(dt,exvCF,imuMeas{i},Pin,biasEstimate,L_ab);
    end
    
    zk=[exvCF(1:6)-gpsMeas{1}(2:7)];
    eye3=eye(3); zer3=zeros(3,3);

    %note: biasState = [de, dx, dv, ba, bg]
    thisEuler=ewxv(1:3);
    RR=calculate_R_from_euler(ewxv(1:3))';
%     %modified Groves formulation
% %     F21=hatmat(RR*imuMeas{numImuMeas}(2:4));
% %     F23=zeros(3,3); F23(3,3) = -2*9.81*Re^2/(Re+ewxv(9))^3;  %effect of deltaz on gravity
%     F21=zeros(3,3); F23=zeros(3,3);
%     Fk=[eye3 zer3 zer3 zer3 RR*dtt
%         zer3 eye3 eye3*dtt RR*dtt^2/2 zer3
%         F21 zer3 F23*dtt RR*dtt zer3
%         zer3 zer3 zer3 eye3 zer3
%         zer3 zer3 zer3 zer3 eye3];
    %modified Brown and Hwang formulation
    Gt=[zer3 dtt*eye3
       RR*dtt^2/2*eye3 zer3
       RR*dtt*eye3 zer3
       eye3 zer3
       zer3 eye3]; %Gammak    
    Fb = hatmat(imuMeas{numImuMeas}(5:7));
    Omegab=hatmat(imuMeas{numImuMeas}(2:4));
    L1b=hatmat(L_ab);
    Fk23 = -RR*Fb;
    Fk33 = -Omegab;
    Ak=[zer3 eye3 zer3 zer3 zer3
        zer3 zer3 Fk23 zer3 RR
        zer3 zer3 Fk33 eye3 zer3
        zeros(6,15)];
    Fk=eye(15)+dtt*Ak;
    zk=[exvCF(4:6)-gpsMeas{1}(5:7)];  %pose meas in local frame
    Rk=.1*eye(3);
    Hk=[eye3 zer3 -RR*L1b zer3 zer3];
%     Rk=.05*eye(6);
%     Hk=[eye3 zer3 zer3 zer3 zer3
%         zer3 eye3 zer3 zer3 zer3];
    biaskp1=Fk*biasState; %biaskp1(6)=biaskp1(6)+9.81*dtt^2/2; biaskp1(9)=biaskp1(9)+9.81*dtt;
    FPF=Fk*Pin*Fk';
    Pkp1=Fk*Pin*Fk'+Gt*Qimu*Gt';
    Sk=Rk+Hk*Pkp1*Hk';
    Kg=Pkp1*Hk'*inv(Sk);
    nuj=zk-Hk*biaskp1;
    biasState=biaskp1+Kg*nuj;
    
end

if numImuMeas>1 %split for initialization
    %ewxvUp=[exvCF(1:3);imuMeas{numImuMeas}(2:4);-biasState(4:6)+gpsMeas{1}(5:7);-biasState(7:9)];
    ewxvUp=[exvCF(1:3);imuMeas{numImuMeas}(2:4);exvCF(4:9)];
    Pout=Pkp1-Kg*Hk*Pkp1;
else
    ewxvUp=ewxv;
    Pout=Pin;
end

end

