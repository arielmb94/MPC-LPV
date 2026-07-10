function mpc = equality_residuals(mpc)

% k = 0
mpc.rp_0(1:mpc.nx) = mpc.B*mpc.u(:,1)-mpc.s(:,1)-mpc.beq_0;

for k = 1:mpc.N-1
    % [A B -I][s u s+]'-beq
    mpc.rp_k(1:mpc.nx,k) = mpc.A*mpc.s(:,k)+mpc.B*mpc.u(:,k+1)-mpc.s(:,k+1)...
                            -mpc.beq_k(:,k);
end

if mpc.has_du
    % k = 0
    mpc.rp_0(mpc.nx+1:mpc.nse) = mpc.u(:,1)-mpc.su(:,1);
    % [I -I][u su+]'
    mpc.rp_k(mpc.nx+1:mpc.nse,:) = mpc.u(:,2:mpc.N)-mpc.su(:,2:mpc.N);
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
for k = 1:mpc.N-1

    mpc.ri_k(:,k) = mpc.Ai_k(:,s_col,k)*mpc.s(:,k)+...
                        mpc.Ai_k(:,u_col,k)*mpc.u(:,k+1)+mpc.g_k(:,k)-...
                        mpc.bi_k(:,k);

    if mpc.has_du
        mpc.ri_k(:,k) = mpc.ri_k(:,k) + mpc.Ai_k(:,su_col,k)*mpc.su(:,k);
    end
    if any(mpc.v_index_k)
        mpc.ri_k(mpc.vi_k,k) = mpc.ri_k(mpc.vi_k,k) - mpc.v_k(:,k);
    end

end

if mpc.has_s_cnstr
    mpc.ri_ter(:) = mpc.Ai_ter(:,s_col)*mpc.s(:,mpc.N)+mpc.g_ter-mpc.v_ter-mpc.bi_ter;
end

end