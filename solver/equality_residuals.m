function mpc = equality_residuals(mpc)

% k = 0
mpc.rp_0(mpc.s_col) = mpc.B*mpc.u(:,1)-mpc.s(:,1)-mpc.beq_0;

for k = 1:mpc.N-1
    % [A B -I][s u s+]'-beq
    mpc.rp_k(mpc.s_col,k) = mpc.A*mpc.s(:,k)+mpc.B*mpc.u(:,k+1)-mpc.s(:,k+1)...
                            -mpc.beq_k(:,k);
end

if mpc.has_du
    % this is the control input delay equality condition: su+ = u_prev
    % [I -I][u su+]' = 0
    % k = 0
    mpc.rp_0(mpc.su_col) = mpc.u(:,1)-mpc.su(:,1);
    % k = 1...N
    mpc.rp_k(mpc.su_col,:) = mpc.u(:,2:mpc.N)-mpc.su(:,2:mpc.N);
end

%%
% k = 0
mpc.ri_0(:) = mpc.Ai_0*mpc.u(:,1)+mpc.g_0-mpc.bi_0;
if mpc.nv_k(1)
    mpc.ri_0(mpc.v_rows_0) = mpc.ri_0(mpc.v_rows_0) - mpc.v_0;
end

for k = 1:mpc.N-1
    % ri = Ai*[s su u]'+ g - v - bi
    mpc.ri_k(:,k) = mpc.Ai_k(:,mpc.s_col,k)*mpc.s(:,k)+...
                        mpc.Ai_k(:,mpc.u_col,k)*mpc.u(:,k+1)+mpc.g_k(:,k)-...
                        mpc.bi_k(:,k);

    if mpc.has_du
        mpc.ri_k(:,k) = mpc.ri_k(:,k) + mpc.Ai_k(:,mpc.su_col,k)*mpc.su(:,k);
    end
    if mpc.nv_k(2)
        mpc.ri_k(mpc.v_rows_k,k) = mpc.ri_k(mpc.v_rows_k,k) - mpc.v_k(:,k);
    end

end

if mpc.ng_k(3)
    % ri = Ai_ter*s_N + g_ter - v_ter - bi_ter
    mpc.ri_ter(:) = mpc.Ai_ter(:,mpc.s_col)*mpc.s(:,mpc.N)+mpc.g_ter-mpc.v_ter-mpc.bi_ter;
end

end