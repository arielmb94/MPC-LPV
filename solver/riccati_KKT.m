function mpc = riccati_KKT(mpc,Q,Q_ter,R_0,R,Y,...
                            ru_0,ru,rs,rs_ter,rp_0,rp)

Q_hat(:,:,mpc.N) = Q_ter;
rs_hat(:,mpc.N) = rs_ter;

for k = mpc.N-1:-1:1

    % \hat{R}_k = R_k + B^T\hat{Q}_{k+1}B
    R_hat(:,:,k) = R(:,:,k) + mpc.B_kkt'*Q_hat(:,:,k+1)*mpc.B_kkt;
    iR(:,:,k) = inv(R_hat(:,:,k));
    % \hat{Y}_k = Y_k + B^T\hat{Q}_{k+1}A
    Y_hat(:,:,k) = Y(:,:,k) + mpc.B_kkt'*Q_hat(:,:,k+1)*mpc.A_kkt;
    % \hat{Q}_k = Q_k + A^T\hat{Q}_{k+1}A - \hat{Y}_k^T\hat{R}_k^{-1}\hat{Y}_k
    Q_hat(:,:,k) = Q(:,:,k) + mpc.A_kkt'*Q_hat(:,:,k+1)*mpc.A_kkt-Y_hat(:,:,k)'*iR(:,:,k)*Y_hat(:,:,k);

    % \hat{r}_{pk} = \hat{r}_{sk+1} + \hat{Q}_{k+1}r_{pk}
    rp_hat(:,k) = rs_hat(:,k+1) + Q_hat(:,:,k+1)*rp(:,k);
    % \hat{r}_{uk} = r_{uk} + B^T\hat{r}_{pk}
    ru_hat(:,k) = ru(:,k) + mpc.B_kkt'*rp_hat(:,k);
    % \hat{r}_{sk} = r_{sk} + A^T\hat{r}_{pk} - \hat{Y}_k^T\hat{R}_k^{-1}\hat{r}_{uk}
    rs_hat(:,k) = rs(:,k) + mpc.A_kkt'*rp_hat(:,k) - Y_hat(:,:,k)'*iR(:,:,k)*ru_hat(:,k);
end

R_hat_0 = R_0+mpc.B_kkt'*Q_hat(:,:,1)*mpc.B_kkt;

rp_hat_0 = rs_hat(:,1) + Q_hat(:,:,1)*rp_0;
ru_hat_0 = ru_0 + mpc.B_kkt'*rp_hat_0;

% u_0 = -\hat{R}_0^{-1}\hat{r}_{u0}
mpc.delta_u(:,1) = -inv(R_hat_0)*ru_hat_0;
% s_1 = Bu_0 + r_{p0}
mpc.delta_se(:,1) = mpc.B_kkt*mpc.delta_u(:,1) + rp_0;
% \mu_k = \hat{r}_{sk+1} + \hat{Q}_{k+1}s_{k+1}
%mu(:,1) = rs_hat(:,1) + Q_hat(:,:,1)*mpc.delta_se(:,1);

for k = 1:mpc.N-1
    ku = k+1; % u mu vectors shifted by 1 stage
    % u_k = -\hat{R}_k^{-1}(\hat{r}_{uk} + \hat{Y}_ks_k)
    mpc.delta_u(:,ku) = -iR(:,:,k)*ru_hat(:,k) - iR(:,:,k)*Y_hat(:,:,k)*mpc.delta_se(:,k);
    % s_{k+1} =  As_k + Bu_k + r_{pk}
    mpc.delta_se(:,k+1) = mpc.A_kkt*mpc.delta_se(:,k) + mpc.B_kkt*mpc.delta_u(:,ku) + rp(:,k);
    % \mu_k = \hat{r}_{sk+1} + \hat{Q}_{k+1}s_{k+1}
    %mu(:,ku) = rs_hat(:,k+1) + Q_hat(:,:,k+1)*mpc.delta_se(:,k+1);
end

end


