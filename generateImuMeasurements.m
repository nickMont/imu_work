function [imuMeas,imuInternalState] = generateImuMeasurements(Rimu,imuInternalState,t_ode,uprate_imu,imuConsts,t_ewxv_ode,ewxv_ode)
%NO LONGER COMPATIBLE
twhole=floor(t_ode(1));
tfrac0=mod(t_ode(1),1);
tfracf=mod(t_ode(2),1);
if tfracf<tfrac0 %handles wraparound
    tfracf=tfracf+1;
end
imutimer=0:1/uprate_imu:2;  %imu has true time twhole+imutimer(timeIndex)
timu_ind=find(imutimer<=tfracf & imutimer>=tfrac0);
imuMeas={}; %cell array of 7x1 measurements, where 1->time, 2:7->meas
cholR=chol(Rimu);

% %set to zero to test noiseless condition. Hardcode rather than flag
% cholR=zeros(9,9);

%ewxvDotMat=zeros(12,length(timu_ind));
ewxvMat=(interp1(t_ewxv_ode,ewxv_ode,twhole+imutimer(timu_ind)))';
%imutimer(timu_ind)
%size(ewxvMat)

for j=1:length(timu_ind)
    %propagate biases
    [imuInternalState,accB,gyroB]=imuPropagateBiasesDiscrete(imuInternalState,1/uprate_imu,imuConsts,cholR*randn(9,1));

    %Sample system by index matching (inaccurate) or interpolation (most accurate)
    % %index matching
    %[~,mindex]=min(abs(t_ewxv_ode-(twhole+imutimer(timu_ind(j)))));
    %ewxvDot=eI_wI_x_v_dot(t_ewxv_ode(mindex),ewxv_ode(mindex,:)');
    %grav=calculate_R_from_euler(ewxv_ode(mindex,1:3))*[0;0;-9.81];

    %interpolation
%    this_ewxv=(interp1(t_ewxv_ode,ewxv_ode,twhole+imutimer(timu_ind(j))))';
    this_ewxv=ewxvMat(:,j);
    ewxvDot=eI_wI_x_v_dot(twhole+imutimer(timu_ind(j)),this_ewxv);  %global frame properties
    RR=calculate_R_from_euler(this_ewxv(1:3));
    grav=[0;0;-9.81];
    
    %combine bias and ewxvDot
    imuMeas{j}=[twhole+imutimer(timu_ind(j));RR*(ewxvDot(1:3))+gyroB;RR*(ewxvDot(10:12)+grav)+accB];
end

end

