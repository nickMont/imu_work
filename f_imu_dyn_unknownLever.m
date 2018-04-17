function xkp1 = f_imu_dyn_unknownLever(dt, xk, RBI0, fB_meas, wB_meas, vk, ...
    tau_accel, tau_gyro,L0)
%xk = [rI; vI; gammak; ba; bg; rIMU];
%vk = [gryo_output; gyro_bias; accel_output; accel_bias];
rI = xk(1:3);
vI = xk(4:6);
gammak = xk(7:9);
ba = xk(10:12);
bg = xk(13:15);
rIMU = xk(16:18)+L0;

RR = euler2dcm(gammak)*RBI0;

wB = wB_meas - bg - vk(1:3);
wBxwBxL = cross(wB,cross(wB,rIMU));
%wBxwBxL = zeros(3,1);
a_est = RR'*(fB_meas - ba - vk(7:9) - wBxwBxL) + [0;0;-9.81];

rIkp1 = rI + vI*dt + a_est*dt^2/2;
vIkp1 = vI + a_est*dt;
gammakp1 = gammak + wB*dt;
bakp1 = (1-dt/tau_accel)*ba + vk(10:12);
bgkp1 = (1-dt/tau_gyro)*bg + vk(4:6);
rIMUkp1 = xk(16:18);

xkp1 = [rIkp1;vIkp1;gammakp1;bakp1;bgkp1;rIMUkp1];

end

