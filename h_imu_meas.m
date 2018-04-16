function zk = h_imu_meas(xk, RBI, vk, rcg2p, rs2p)

zk=zeros(6,1);
RR = euler2dcm(xk(7:9))*RBI;
zk(1:3) = xk(1:3)+RR'*rcg2p+vk(1:3);
zk(4:6) = RR'*unit3(rs2p+vk(4:6));

end

