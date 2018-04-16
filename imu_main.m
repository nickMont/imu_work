clear;clc;
% Generates and samples IMU data

% GPS/INS sampling frequencies
fratio=4; %Ratio of imu to gps sampling rates. INTEGER ONLY
fgps=20;
fimu=fratio*fgps;

% GPS parameters
Rgps(:,:,1)=.02*eye(3);
Rgps(:,:,2)=Rgps(:,:,1);
cholR1=chol(Rgps(:,:,1));
cholR2=chol(Rgps(:,:,2));
antenna2_loc = [-1;0;0];

% INS noise parameters

% Before subsampling
thist0=0:1/fgps:20;
N0=length(thist0);
xhist_init=zeros(401,3);  %position hist in inertial frame
xhist_init(:,1)=sin(thist0);
ehist_init=zeros(401,3);  %attitude hist for 312 convention

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
Rimu=diag([1e-12*ones(3,1); 1000*(1/9.81*10^-6)^2*ones(3,1); 10^-8*(pi/180)^2*ones(3,1)])*fimu;
imuConsts=[100;100;100;10;10;10;10;10;10]; %imu decorrelation times, see imuPropagateBiasesDiscrete
imuInternalState=zeros(12,1);
statehist=[xhist vhist ahist ehist omegaBhist alphaBhist];
Lab_true=[0;0;0];
for ij=2:N0
    z1 = xhist_init(ij,:)' + cholR1*randn(3,1);
    z2 = xhist_init(ij,:)' + cholR2*randn(3,1)+...
        calculate_R_from_euler(ehist_init(ij,:))*antenna2_loc;
    gpsMeas = {[thist0(ij);z1;z2]};
    gpsMeasStore{ij-1}=gpsMeas;
    
%     for ik=1:5
%         iI=(ij-2)*fratio+ik;
%         zImu(1)=thist(iI);
%         imuMeas{ik}=zIImu;
%     end
    t_ode=[thist((ij-2)*fratio+1) : 1/fimu : thist((ij-2)*fratio+fratio+1)];
    %t_ode=[thist((ij-2)*fratio+1) thist((ij-2)*fratio+fratio+1)];
    [imuMeas,imuInternalState]=generateImuMeasurementsWithLever(Rimu,...
        imuInternalState, t_ode, fimu, imuConsts, thist, statehist, Lab_true);
    %imuMeas{1}
    imuMeasStore{ij-1}=imuMeas;
end
%NOTE: gpsMeasStore{1} corresponds to thist0(2)











