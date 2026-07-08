function [u,s,mu] = riccati_KKT(mpc,Q,R,Y,ru,rs,rp)

Q_hat(:,:,mpc.N) = Q(:,:,mpc.N);
rs_hat(:,mpc.N) = rs(:,mpc.N);

for k = mpc.N-1:-1:1

    % \hat{R}_k = R_k + B^T\hat{Q}_{k+1}B
    R_hat(:,:,k) = R(:,:,k+1) + mpc.B_kkt'*Q_hat(:,:,k+1)*mpc.B_kkt;
    iR(:,:,k) = inv(R_hat(:,:,k));
    % \hat{Y}_k = Y_k + B^T\hat{Q}_{k+1}A
    Y_hat(:,:,k) = Y(:,:,k) + mpc.B_kkt'*Q_hat(:,:,k+1)*mpc.A_kkt;
    % \hat{Q}_k = Q_k + A^T\hat{Q}_{k+1}A - \hat{Y}_k^T\hat{R}_k^{-1}\hat{Y}_k
    Q_hat(:,:,k) = Q(:,:,k) + mpc.A_kkt'*Q_hat(:,:,k+1)*mpc.A_kkt-Y_hat(:,:,k)'*iR(:,:,k)*Y_hat(:,:,k);

    % \hat{r}_{pk} = \hat{r}_{sk+1} + \hat{Q}_{k+1}r_{pk}
    rp_hat(:,k) = rs_hat(:,k+1) + Q_hat(:,:,k+1)*rp(:,k+1);
    % \hat{r}_{uk} = r_{uk} + B^T\hat{r}_{pk}
    ru_hat(:,k) = ru(:,k+1) + mpc.B_kkt'*rp_hat(:,k);
    % \hat{r}_{sk} = r_{sk} + A^T\hat{r}_{pk} - \hat{Y}_k^T\hat{R}_k^{-1}\hat{r}_{uk}
    rs_hat(:,k) = rs(:,k) + mpc.A_kkt'*rp_hat(:,k) - Y_hat(:,:,k)'*iR(:,:,k)*ru_hat(:,k);
end

R0 = R(:,:,1)+mpc.B_kkt'*Q_hat(:,:,1)*mpc.B_kkt;

rp_hat_0 = rs_hat(:,1) + Q_hat(:,:,1)*rp(:,1);
ru_hat_0 = ru(:,1) + mpc.B_kkt'*rp_hat_0;

% u_0 = -\hat{R}_0^{-1}\hat{r}_{u0}
u(:,1) = -inv(R0)*ru_hat_0;
% s_1 = Bu_0 + r_{p0}
s(:,1) = mpc.B_kkt*u(:,1) + rp(:,1);
% \mu_k = \hat{r}_{sk+1} + \hat{Q}_{k+1}s_{k+1}
mu(:,1) = rs_hat(:,1) + Q_hat(:,:,1)*s(:,1);

for k = 2:mpc.N
    % u_k = -\hat{R}_k^{-1}(\hat{r}_{uk} + \hat{Y}_ks_k)
    u(:,k) = -iR(:,:,k-1)*ru_hat(:,k-1) - iR(:,:,k-1)*Y_hat(:,:,k-1)*s(:,k-1);
    % s_{k+1} =  As_k + Bu_k + r_{pk}
    s(:,k) = mpc.A_kkt*s(:,k-1) + mpc.B_kkt*u(:,k) + rp(:,k);
    % \mu_k = \hat{r}_{sk+1} + \hat{Q}_{k+1}s_{k+1}
    mu(:,k) = rs_hat(:,k) + Q_hat(:,:,k)*s(:,k);
end

end


