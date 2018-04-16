function [xbarkp1,Pbarkp1] = ukfPropagate(dt,xk,Pk,Qk, RBI0,fB,wB,Limu0,alphaA,alphaG)

nx = length(xk);
[~,nv] = size(Qk);
alpha = 1e-3;
beta = 2;
kappa = 0; 
lambda_p = alpha^2*(kappa + nx + nv) - nx - nv;
c_p = sqrt(nx+nv+lambda_p);
w_mean_center = lambda_p/(nx + nv + lambda_p);
w_mean_reg = 1/(2*(nx + nv + lambda_p));
w_cov_center = w_mean_center + 1 - alpha^2 + beta;
w_cov_reg = w_mean_reg;

xhat_aug = [xk; zeros(nv,1)];
P_aug = blkdiag(Pk,Qk);
Sk_aug = chol(P_aug)';

% Propagate regression points
sp0 = xHatAugk;
xpMat = zeros(nx,1+2*(nx+nv));
xpMat(:,1) = f_imu_dyn(dt,xhat_aug(1:nx),RBI0,fB,wB,xhat_aug(nx+1:end),alphaA,alphaG,Limu0);
xbarkp1 = xpMat(:,1)*w_mean_center;
sgn = 1;
for ij=1:2*(nx+nv)
    colno = mod(ij,nx+nv)+1;
    if(ij > (nx + nv))
        sgn = -1;
    end
    xaug_ij = sp0 + sgn*c_p*Sk_aug(:,colno);
    xpMat(:,ij+1) = f_imu_dyn(dt,xaug_ij(1:nx),RBI0,fB,wB,xaug_ij(nx+1:end),alphaA,alphaG,Limu0);
    xbarkp1 = xbarkp1 + w_mean_reg*xpMat(:,ij+1);
end

%Form covariance
Pbarkp1 = w_cov_center*(xpMat(:,1) - xbarkp1)*(xpMat(:,1) - xbarkp1)';
for ij=1:2*(nx+nv)
    Pbarkp1 = Pbarkp1 + w_cov_reg*(xpMat(:,ij+1) - xbarkp1)*(xpMat(:,ij+1) - xbarkp1)';
end

end

