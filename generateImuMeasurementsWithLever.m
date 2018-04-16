function [imuMeas,imuInternalState] = generateImuMeasurementsWithLever(Rimu,imuInternalState,t_ode,fimu,imuConsts,t_samp,state_samp,L_ab)
%Takes in state as [x v a e wB aB]

twhole=floor(t_ode(1));
tfrac0=mod(t_ode(1),1);
tfracf=mod(t_ode(end),1);
if tfracf<tfrac0 %handles wraparound
    tfracf=tfracf+1;
end
imutimer=0:1/fimu:2;  %imu has true time twhole+imutimer(timeIndex)
timu_ind=find(imutimer>=(tfrac0-1e-7) & imutimer<=(tfracf+1e-7));
imuMeas={}; %cell array of 7x1 measurements, where 1->time, 2:7->meas
cholR=chol(Rimu);

stateMat=(interp1(t_samp,state_samp,twhole+imutimer(timu_ind)))';
%imutimer(timu_ind)

for ij=1:length(timu_ind)
    %propagate biases
    [imuInternalState,accB,gyroB]=imuPropagateBiasesDiscrete(imuInternalState,1/fimu,imuConsts,cholR*randn(9,1));

    %Sample system by index matching (inaccurate) or interpolation (most accurate)
    % %index matching
    %[~,mindex]=min(abs(t_ewxv_ode-(twhole+imutimer(timu_ind(j)))));
    %ewxvDot=eI_wI_x_v_dot(t_ewxv_ode(mindex),ewxv_ode(mindex,:)');
    %grav=calculate_R_from_euler(ewxv_ode(mindex,1:3))*[0;0;-9.81];

    %interpolation
%    this_ewxv=(interp1(t_ewxv_ode,ewxv_ode,twhole+imutimer(timu_ind(j))))';
    this_ewxv=stateMat(:,ij);
    RR=calculate_R_from_euler(this_ewxv(10:12));
    wB=this_ewxv(13:15);
    grav=[0;0;9.81];
   
    outputNoise = zeros(6,1);
    
    %combine bias and ewxvDot
    imuMeas{ij}=[twhole+imutimer(timu_ind(ij));
        wB+gyroB+outputNoise(1:3);
        RR*(this_ewxv(7:9) + grav + accB) + (cross(this_ewxv(16:18),L_ab)+cross(wB,cross(wB,L_ab)))+outputNoise(4:6)];
end

end

