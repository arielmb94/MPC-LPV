function mpc = equality_residuals(mpc)

% k = 0
mpc.rp(1:mpc.nx,1) = mpc.B*mpc.u(:,1)-mpc.s(:,2)-mpc.beq(:,1);

for k = 2:mpc.N
    mpc.rp(1:mpc.nx,k) = mpc.A*mpc.s(:,k)+mpc.B*mpc.u(:,k)-mpc.s(:,k+1)...
                            -mpc.beq(:,k);
end

if mpc.has_du
    for k = 1:mpc.N
        mpc.rp(mpc.nx+1:mpc.nse,k) = mpc.u(:,k)-mpc.su(:,k+1);
    end
end

%%
% k = 0
mpc.ri_0(:) = mpc.Ai_0*mpc.u(:,1)+mpc.g_0-mpc.bi_0;
if any(mpc.v_index_0)
    mpc.ri_0(mpc.vi_0) = mpc.ri_0(mpc.vi_0) - mpc.v_0;
end

s_col = mpc.inq_s_col;
su_col = mpc.inq_su_col;
u_col = mpc.inq_u_col;
for k = 2:mpc.N
    ki = k-1; % inequalities from k = 1,...,N-1 are shifted with respect variables
    mpc.ri_k(:,ki) = mpc.Ai_k(:,s_col,ki)*mpc.s(:,k)+...
                        mpc.Ai_k(:,u_col,ki)*mpc.u(:,k)+mpc.g_k(:,ki)-...
                        mpc.bi_k(:,ki);
    if mpc.has_du
        mpc.ri_k(:,ki) = mpc.ri_k(:,ki) + mpc.Ai_k(:,su_col,ki)*mpc.su(:,k);
    end
    if any(mpc.v_index_k)
        mpc.ri_k(mpc.vi_k,ki) = mpc.ri_k(mpc.vi_k,ki) - mpc.v_k(:,ki);
    end

end

if mpc.has_s_cnstr
    mpc.ri_ter(:) = mpc.Ai_ter(:,s_col)*mpc.s(:,mpc.N+1)+mpc.g_ter-mpc.v_ter-mpc.bi_ter;
end

end