function [zbar,Pxz,Pzz] = ukfMeasure(xbar,Pbar,Rk, RBI0,Lcg2p,Ls2p)

nx = length(xbar);
[~,nw] = size(Rk);
nz=6;
alpha = 1e-3;
beta = 2;
kappa = 0; 
lambda_p = alpha^2*(kappa + nx + nw) - nx - nw;
c_p = sqrt(nx+nw+lambda_p);
w_mean_center = lambda_p/(nx + nw + lambda_p);
w_mean_reg = 1/(2*(nx + nw + lambda_p));
w_cov_center = w_mean_center + 1 - alpha^2 + beta;
w_cov_reg = w_mean_reg;

xhat_aug = [xbar; zeros(nw,1)];
P_aug = blkdiag(Pbar,Rk);
Sk_aug = chol(P_aug)';

% Propagate regression points
zpMat = zeros(nz,1+2*(nx+nw));
zpMat(:,1) = h_imu_meas(xhat_aug(1:nx),RBI0,xhat_aug(nx+1:end),Lcg2p,Ls2p);
zbar = zpMat(:,1)*w_mean_center;
sgn = 1;
xpMat=zeros(nx,1+2*(nx+nw));
xpMat(:,1) = xhat_aug(1:nx);
for ij=1:2*(nx+nw)
    colno = mod(ij,nx+nw)+1;
    if(ij > (nx + nw))
        sgn = -1;
    end
    xaug_ij = xhat_aug + sgn*c_p*Sk_aug(:,colno);
    xpMat(:,ij+1) = xaug_ij(1:nx);
    zpMat(:,ij+1) = h_imu_meas(xaug_ij(1:nx),RBI0,xaug_ij(nx+1:end),Lcg2p,Ls2p);
    zbar = zbar + w_mean_reg*zpMat(:,ij+1);
end

%Form covariance
Pzz = w_cov_center*(zpMat(:,1) - zbar)*(zpMat(:,1) - zbar)';
Pxz = w_cov_center*(xpMat(:,1)-xbar)*(zpMat(:,1)-zbar)';
for ij=1:2*(nx+nw)
    Pzz = Pzz + w_cov_reg*(zpMat(:,ij+1) - zbar)*(zpMat(:,ij+1) - zbar)';
    Pxz = Pxz + w_cov_reg*(xpMat(:,ij+1)-xbar)*(zpMat(:,ij+1)-zbar)';
end

end

