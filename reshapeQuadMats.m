clear;clc;
loadline=44;
initstate=[];

rcg2p=[0;0;0];

load dat999.mat
data=squeeze(trainDat(loadline,:,1,1));
ewxv0=data(1:12);
x0=ewxv0(7:9);
v0=ewxv0(10:12);
RBI0=euler2dcm(ewxv0(1:3));
gpsMeasStore={};
imuMeasStore={};

for ijk=1:400
    gpsMeas={data(12+(ijk-1)*4+1:12+4*ijk)'};
    gpsMeasStore{ijk}=gpsMeas;
end
imuData=data(12+4*400+1:end)';
times=imuData(1:7:end);

tlast=0;
for ijk=1:400
    gpsTime=gpsMeasStore{ijk}{1}(1);
    imuInd=find((times(:)>(tlast-1e-9)) & (times(:)<(gpsTime+1e-9)));
    for j=1:length(imuInd)
        tt=imuInd(j);
        imuMeas{j}=imuData((tt-1)*7+1:tt*7);
    end
    imuMeasStore{ijk}=imuMeas;
    
end

% GPS/INS sampling frequencies
fratio=4; %Ratio of imu to gps sampling rates. INTEGER ONLY
fgps=20;
fimu=fratio*fgps;

% GPS parameters
Rgps(:,:,1)=.0002*eye(3);
Rgps(:,:,2)=Rgps(:,:,1);
cholR1=chol(Rgps(:,:,1));
cholR2=chol(Rgps(:,:,2));
%cholR1=zeros(3,3);
%cholR2=zeros(3,3);
rs2p = [1;0;0];
rcg2p=[0;0;0];

%imu consts
imuConsts=[100;100;100;100;100;100];
Rimu.bias = diag([1000*(1/9.81*10^-6)^2*ones(3,1); 10^-8*(pi/180)^2*ones(3,1)])*fimu;
Rimu.output = diag([(9.81/1e5/sqrt(10))^2*ones(3,1); 25*(pi/180/100)^2*ones(3,1)])*fimu;


% Match convention in f_imu_dyn
RimuStack = blkdiag(Rimu.output(4:6,4:6),Rimu.bias(4:6,4:6),Rimu.output(1:3,1:3),Rimu.bias(1:3,1:3));
systemParams.Rimu = RimuStack;
systemParams.tA = imuConsts(1);
systemParams.tG = imuConsts(4);
systemParams.Ls2p = rs2p;
systemParams.Lcg2p = rcg2p;

RBI = euler2dcm(ewxv0(1:3));
limu0=[0;0;0];
%x0=xhist(1,:)';
%v0=vhist(1,:)';
g0=zeros(3,1);
ba0=zeros(3,1);
bg0=zeros(3,1);

nmax = length(gpsMeasStore);

% UKF with state augmentation
state0=[x0';v0';g0;ba0;bg0;limu0];
P0=diag([.1*ones(3,1); .01*ones(3,1); .01*ones(3,1); ...
    .01*ones(3,1); .01*ones(3,1); 0.05*ones(3,1)]);
statestore=zeros(18,nmax);

% %15
% state0=[x0;v0;g0;ba0;bg0];
% P0=diag([.1*ones(3,1); .01*ones(3,1); .01*ones(3,1); ...
%     .01*ones(3,1); .01*ones(3,1)]);
% limu0=[-1;0;0];
% statestore=zeros(15,nmax);


% state aug UKF
for ij=1:nmax
    [state,Pk,RBI,Limu]=runUKF(imuMeasStore{ij},gpsMeasStore{ij},state0,RBI,P0,systemParams);
    %[state,Pk,RBI]=runUKF15(imuMeasStore{ij},gpsMeasStore{ij},state0,RBI,P0,systemParams,limu0);
    Pk=Pk+1e-10*eye(length(Pk));
    LL=Limu
    statestore(:,ij)=state;
    %dxv=state(1:6)-[xhist_init(ij+1,:)';vhist(4*(ij-1)+1,:)']
    %interestingStates=[state(1:6);state(10:15)]
    RBI
end