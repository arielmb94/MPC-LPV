function mpc = recover_slacks(mpc,delta_u,delta_se)
%% \mu_i = (-S)^{-1}(-\hat{r}_i-A_i\Delta x) = S^{-1}(\hat{r}_i+A_i\Delta x)

%k = 0
mpc.mu_i_0(:) = mpc.iS_ri_hat_0+mpc.iS_Ai_0*delta_u(:,1);

for k = 1:mpc.N-1
    mpc.mu_i_k(:,k) = mpc.iS_ri_hat_k(:,k) +...
                      mpc.iS_Ai_k(:,mpc.se_col,k)*delta_se(:,k) +...
                      mpc.iS_Ai_k(:,mpc.u_col,k)*delta_u(:,k+1);
end
%k = N
if mpc.ng_k(3)
    mpc.mu_i_ter(:) = mpc.iS_ri_hat_ter+mpc.iS_Ai_ter*delta_se(:,mpc.N);
end

%% \Delta g  = (g^2)(-r_g-\mu_i) 
%  \Delta v = (v^2)(-r_v+ \mu_i) 

mpc.delta_g_0(:) = mpc.g_0 - (mpc.g_0.^2).*mpc.mu_i_0;
if mpc.nv_k(1)
    mpc.delta_v_0(:) = -mpc.rv_0+(mpc.v_0.^2).*mpc.mu_i_0(mpc.v_rows_0);
end

mpc.delta_g_k = mpc.g_k - (mpc.g_k.^2).*mpc.mu_i_k;
if mpc.nv_k(2)
    mpc.delta_v_k(:) = -mpc.rv_k+(mpc.v_k.^2).*mpc.mu_i_k(mpc.v_rows_k,:);
end

if mpc.ng_k(3)
    mpc.delta_g_ter = mpc.g_ter - (mpc.g_ter.^2).*mpc.mu_i_ter;
    mpc.delta_v_ter = -mpc.rv_ter +(mpc.v_ter.^2).*mpc.mu_i_ter;
end

end