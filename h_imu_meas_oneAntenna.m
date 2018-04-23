function zk = h_imu_meas_oneAntenna(xk, RBI, vk, rcg2p)

zk=zeros(3,1);
RR = euler2dcm(xk(7:9))*RBI;
zk(1:3) = xk(1:3)+RR'*rcg2p+vk(1:3);

end

