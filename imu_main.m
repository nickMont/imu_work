clear all;clc;
% Generates and samples IMU data

rng(10)

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

% INS noise parameters

% Before subsampling
thist0=0:1/fgps:100;
N0=length(thist0);
%N0=3;
xhist_init=zeros(N0,3);  %position hist in inertial frame
xhist_init(:,1)=sin(thist0);
ehist_init=zeros(N0,3);  %attitude hist for 312 convention
ehist_init(:,1)=cos(0.5*thist0);
ehist_init(:,3)=0.01*sin(0.1*thist0);

% New sampled time
thist=(thist0(1):1/fimu:thist0(end))';

% Sampled values
N=length(thist);
xhist=zeros(N,3);
vhist=zeros(N,3);
ahist=zeros(N,3);
ehist=zeros(N,3);
edothist=zeros(N,3);
omegaBhist=zeros(N,3);
alphaBhist=zeros(N,3);

% Generate truth data history
dt=thist(2)-thist(1);
window=8;  %size of window for moving average
order=4;   %approximation order for moving average
% Generate x,v,a,e for subsampling
for ij=1:3
    s=spline(thist0,xhist_init(:,ij),thist);
    xhist(:,ij)=s;
    v=movingslope(s,window,order,dt);
    vhist(:,ij)=v;
    ahist(:,ij)=movingslope(v,window,order,dt);
    
    s=spline(thist0,ehist_init(:,ij),thist);
    ehist(:,ij)=s;
    edothist(:,ij)=movingslope(s,window,order,dt);
end
% Rotate omegaI to omegaB
for ij=1:N
    omegaBhist(ij,:)=(R_pqr(ehist(ij,:))*edothist(ij,:)')';
end
% Inerpolate whist to get alphahist
for ij=1:3
    alphaBhist(:,ij)=movingslope(omegaBhist(:,ij),window,order,dt);
end

% Generate GPS and IMU data
gpsMeasStore={};
imuMeasStore={};
zImu=zeros(7,1);

% %Two-state bias model
% Rimu=diag([1e-12*ones(3,1); 1000*(1/9.81*10^-6)^2*ones(3,1); 10^-8*(pi/180)^2*ones(3,1)])*fimu;
% imuConsts=[100;100;100;10;10;10;10;10;10]; %imu decorrelation times, see imuPropagateBiasesDiscrete
% imuInternalState=zeros(12,1);

% One-state bias model with output noise
imuConsts=[100;100;100;100;100;100];
Rimu.bias = diag([1000*(1/9.81*10^-6)^2*ones(3,1); 10^-8*(pi/180)^2*ones(3,1)])*fimu;
Rimu.output = diag([(9.81/1e5/sqrt(10))^2*ones(3,1); 25*(pi/180/100)^2*ones(3,1)])*fimu;
%imuInternalState=zeros(6,1);
imuInternalState=(2*eye(6)-fimu*inv(diag(imuConsts)))*chol(Rimu.bias)*randn(6,1); %approximate SS from expm
statehist=[xhist vhist ahist ehist omegaBhist alphaBhist];
Lab_true=[-1;0;0];
tLast=-1;
testmat=[];
for ij=2:N0
    z1 = xhist_init(ij,:)' + cholR1*randn(3,1) +...
        euler2dcm(ehist_init(ij,:))*rcg2p;
    z2 = euler2dcm(ehist_init(ij,:))*rs2p + cholR2*randn(3,1);
    gpsMeas = {[thist0(ij);z1;z2]};
    gpsMeasStore{ij-1}=gpsMeas;
    
    % %Specify times for measurement       
    %t_ode=[thist((ij-2)*fratio+1) : 1/fimu : thist((ij-2)*fratio+fratio+1)];
    % Specify start and end times and let generateImuMeasurements sort it
    t_ode=[thist((ij-2)*fratio+1) thist((ij-2)*fratio+fratio+1)];
    
    
    % Grab imu measurements
    [imuMeas,imuInternalState]=generateImuMeasurementsWithLever(Rimu,...
        imuInternalState, t_ode, fimu, imuConsts, thist, statehist, Lab_true, tLast);
    tLast = t_ode(end);  %tLast is used internally to prevent biases from being double-propagated
    imuMeasStore{ij-1}=imuMeas;
    for ik=1:length(imuMeas)-1
        testmat=[testmat imuMeas{ik}(5)];
        
    end
end
%NOTE: gpsMeasStore{1} corresponds to thist0(2)

% Match convention in f_imu_dyn
RimuStack = blkdiag(Rimu.output(4:6,4:6),Rimu.bias(4:6,4:6),Rimu.output(1:3,1:3),Rimu.bias(1:3,1:3));
systemParams.Rimu = RimuStack;
systemParams.tA = imuConsts(1);
systemParams.tG = imuConsts(4);
systemParams.Ls2p = rs2p;
systemParams.Lcg2p = rcg2p;

RBI = euler2dcm(ehist(1,:));
limu0=[0;0;0];
x0=xhist(1,:)';
v0=vhist(1,:)';
g0=zeros(3,1);
ba0=zeros(3,1);
bg0=zeros(3,1);
state0=[x0;v0;g0;ba0;bg0;limu0];
P0=diag([.1*ones(3,1); .01*ones(3,1); .01*ones(3,1); ...
    .01*ones(3,1); .01*ones(3,1); 0.1*ones(3,1)]);

nmax = length(gpsMeasStore);
%nmax=2;
statestore=zeros(18,nmax);
for ij=1:nmax
    [state,Pk,RBI,Limu]=runUKF(imuMeasStore{ij},gpsMeasStore{ij},state0,RBI,P0,systemParams);
    Pk=Pk+1e-10*eye(length(Pk));
    LL=Limu
    statestore(:,ij)=state;
    dxv=state(1:6)-[xhist_init(ij+1,:)';vhist(4*(ij-1)+1,:)']
    %interestingStates=[state(1:6);state(10:15)]
    RBI
end




