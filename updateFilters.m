function [ewxvUp,biasState,Pout] = updateFilters(dtt,ewxv,Pin,imuMeas,gpsMeas,biasState)
%get all time stamps. GPS is assumed to be the last measurement in terms of
%time
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
        [exvCF,Pin]=complimentaryFilter(dt,exvCF,imuMeas{i},Pin,biasEstimate);
    end
    
    zk=[exvCF(1:6)-gpsMeas{1}(2:7)];
    eye3=eye(3); zer3=zeros(3,3);
    Gt=[zer3 dtt*eye3
       dtt^2/2*eye3 zer3
       dtt*eye3 zer3
       eye3 zer3
       zer3 eye3]; %Gammak
    %note: biasState = [de, dx, dv, ba, bg]
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
    Fk13=dtt/Re*[1 0 0; 0 1 0; 0 0 0];
    Fk31=dtt*-9.81*[1 0 0; 0 1 0; 0 0 0];
    Fk=[eye3 zer3 Fk13 zer3 dtt*RR
        eye3 eye3*dtt zer3 dtt^2/2*RR zer3
        Fk31 zer3 eye3 dtt*RR zer3
        zer3 zer3 zer3 eye3 zer3
        zer3 zer3 zer3 zer3 eye3];
    Rk=.05*eye(6);
    Hk=[eye3 zer3 zer3 zer3 zer3
        zer3 eye3 zer3 zer3 zer3];
    biaskp1=Fk*biasState; %biaskp1(6)=biaskp1(6)+9.81*dtt^2/2; biaskp1(9)=biaskp1(9)+9.81*dtt;
    Pkp1=Fk*Pin*Fk'+Gt*Qimu*Gt';
    Sk=Rk+Hk*Pkp1*Hk';
    Kg=Pkp1*Hk'*inv(Sk);
    biasState=biaskp1+Kg*(zk-Hk*biaskp1);
    
end

if numImuMeas>1 %split for initialization
    ewxvUp=[exvCF(1:3);imuMeas{numImuMeas}(2:4);exvCF(4:9)];
    Pout=Pkp1-Kg*Hk*Pkp1;
else
    ewxvUp=ewxv;
    Pout=Pin;
end

end

