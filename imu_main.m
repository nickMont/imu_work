clear all;clc;
% Generates and samples IMU data

rng(87)

runStateAug=0;
runMMAE=1;
%discSize=0.75; %discretization size for MMAE

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
thist0=0:1/fgps:20;
N0=length(thist0);
%N0=3;
xhist_init=zeros(N0,3);  %position hist in inertial frame
xhist_init(:,1)=sin(thist0);
%xhist_init(:,2)=cos(0.2*thist0);
%xhist_init(:,3)=0.1*thist0;
ehist_init=zeros(N0,3);  %attitude hist for 312 convention
ehist_init(:,1)=cos(0.5*thist0);
ehist_init(:,2)=0.01*cos(0.1*thist0);
ehist_init(:,3)=0.01*sin(0.2*thist0);

Lab_true=[1;0.5;0.2];

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
imu0=imuInternalState;
statehist=[xhist vhist ahist ehist omegaBhist alphaBhist];
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
nmax = length(gpsMeasStore);

datset=[];
datset=[datset statehist(1,:)];
for ij=1:400
    datset=[datset gpsMeasStore{ij}{1}'];
end
for ij=1:400
    n=length(imuMeasStore{ij});
    for iL=1:n
        datset=[datset imuMeasStore{ij}{iL}'];
    end
end
levers=Lab_true';

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

% UKF with state augmentation
state0=[x0;v0;g0;ba0;bg0;limu0];
P0=diag([.1*ones(3,1); .01*ones(3,1); .01*ones(3,1); ...
    .01*ones(3,1); .01*ones(3,1); 0.05*ones(3,1)]);
statestore=zeros(18,nmax);

% %15
% state0=[x0;v0;g0;ba0;bg0];
% P0=diag([.1*ones(3,1); .01*ones(3,1); .01*ones(3,1); ...
%     .01*ones(3,1); .01*ones(3,1)]);
% limu0=[-1;0;0];
% statestore=zeros(15,nmax);

if runStateAug==1
% state aug UKF
for ij=1:nmax
    [state,Pk,RBI,Limu]=runUKF(imuMeasStore{ij},gpsMeasStore{ij},state0,RBI,P0,systemParams);
    %[state,Pk,RBI]=runUKF15(imuMeasStore{ij},gpsMeasStore{ij},state0,RBI,P0,systemParams,limu0);
    Pk=Pk+1e-10*eye(length(Pk));
    LL=Limu;
    statestore(:,ij)=state;
    dxv=state(1:6)-[xhist_init(ij+1,:)';vhist(4*(ij-1)+1,:)'];
    %interestingStates=[state(1:6);state(10:15)]
    RBI;
end

figure(1);clf;
pli=1:1:length(statestore);
hold on
plot(thist0(1+pli),statestore(16,pli)-Lab_true(1),':r')
hold on
plot(thist0(1+pli),statestore(17,pli)-Lab_true(2),'--b')
hold on
plot(thist0(1+pli),statestore(18,pli)-Lab_true(3),'-g')
hold on
legend('\delta l_x','\delta l_y','\delta l_z')
axis([0 20 -1.5 1.5])
xlabel('Time (s)')
ylabel('Error (m)')
grid on
figset
end

if runMMAE==1
% MMAE
ba0=imu0(1:3);
bg0=imu0(4:6);
mu_min=1e-8;
%leverSet=permn(-1:discSize:1,3);
%leverSet=[-2 -2 -2;1 .5 .2];
%numLevers=length(leverSet);
leverSet=(Lab_true*( 0.7:.05:1.3 ))';
[numLevers,~]=size(leverSet);
stateSet=zeros(15,numLevers);
PkSet=zeros(15,15,numLevers);
lambdaTemp=zeros(numLevers,1);
state15=[x0;v0;g0;ba0;bg0];
P15=diag([.1*ones(3,1); .01*ones(3,1); .01*ones(3,1); ...
    .0001*ones(3,1); .0001*ones(3,1)]);
for ijk=1:numLevers
    PkSet(:,:,ijk)=P15;
    stateSet(:,ijk)=state15;
end
mukhist=zeros(numLevers,nmax);
muPrev=ones(numLevers,1)*1/numLevers;
leverEst=zeros(3,nmax);
numActiveLevers=numLevers;
activeLeverSet=1:1:numLevers;
for ij=1:nmax
    tic
    
    %NOTE: Model transition probability is zero
    for ijk=1:numLevers
        [state,Pk,RBI,Sk_notfixed,nuj]=runUKF15(imuMeasStore{ij},gpsMeasStore{ij},stateSet(:,ijk),RBI,PkSet(:,:,ijk),systemParams,leverSet(ijk,:)');
        stateSet(:,ijk)=state;
        %tuning increased gain
        Sk_notfixed = Sk_notfixed+.1*eye(6);
        PkSet(:,:,ijk)=Pk;
        Sk = (Sk_notfixed + Sk_notfixed.')/2; %fix numerical error caused by rounding in inv()
        normpdf_Eval = mvnpdf(nuj,zeros(6,1),Sk);
        lambdaTemp(ijk)=normpdf_Eval;
    end
    %merge to get new mus
    muStackTemp=zeros(numLevers,1);
    for ijk=1:numLevers
       muStackTemp(ijk)=lambdaTemp(ijk)*muPrev(ijk)/dot(lambdaTemp,muPrev); 
    end
    for ijk=1:numLevers
        if muStackTemp(ijk)<=mu_min
            muStackTemp(ijk)=mu_min;
        end
    end
    muTest=muStackTemp/sum(muStackTemp);
    if max(max(isnan(muTest)))>0.1
        error('isnan')
    end
    muPrev=muTest;
     
    mukhist(:,ij)=muPrev;
    L=zeros(3,1);
    for ijk=1:numLevers
        L=L+muPrev(ijk)*leverSet(ijk,:)';
    end
    leverEst(:,ij)=L;
    
    t=toc;
    pctComplete=ij/nmax*100
end

figure(1);clf;
pli=1:1:length(leverEst);
hold on
plot(thist0(1+pli),leverEst(1,pli)-Lab_true(1),':r')
hold on
plot(thist0(1+pli),leverEst(2,pli)-Lab_true(2),'--b')
hold on
plot(thist0(1+pli),leverEst(3,pli)-Lab_true(3),'-g')
hold on
legend('\delta l_x','\delta l_y','\delta l_z')
axis([0 20 -1.5 1.5])
xlabel('Time (s)')
ylabel('Error (m)')
grid on
figset
end


% figure(2);clf;
% hold on
% plot(thist0,xhist_init(:,1),':r')
% hold on
% plot(thist0,xhist_init(:,2),'--b')
% hold on
% plot(thist0,xhist_init(:,3),'-g')
% hold on
% xlabel('Time (s)')
% ylabel('Position (m)')
% legend('x','y','z')
% figset
% % 
% figure(3);clf;
% hold on
% plot(thist0,ehist_init(:,1),':r')
% hold on
% plot(thist0,5*ehist_init(:,2),'--b')
% hold on
% plot(thist0,5*ehist_init(:,3),'-g')
% hold on
% xlabel('Time (s)')
% ylabel('Attitude (radians)')
% legend('Roll','5x Pitch','5x Yaw')
% figset
% 
% figure(4);clf;
% hold on
% plot(thist,vhist(:,1),':r')
% hold on
% plot(thist,vhist(:,2),'--b')
% hold on
% plot(thist,vhist(:,3),'-g')
% hold on
% xlabel('Time (s)')
% ylabel('Velocity (m/s)')
% legend('v_x','v_y','v_z')
% axis([0 20 -1 1])
% figset
% % 
% figure(5);clf;
% hold on
% plot(thist,omegaBhist(:,1),':r')
% hold on
% plot(thist,5*omegaBhist(:,2),'--b')
% hold on
% plot(thist,5*omegaBhist(:,3),'-g')
% hold on
% xlabel('Time (s)')
% ylabel('Attitude (rad/s)')
% legend('\omega_{roll}','5x \omega_{pitch}','5x \omega_{yaw}')
% figset
