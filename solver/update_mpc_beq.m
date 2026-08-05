function mpc = update_mpc_beq(mpc,s_prev,u_prev)
%% dynamcis
% k = 0
mpc.beq_0(:) = -mpc.A(:,:,1)*s_prev;

if mpc.dyn_use_d
    mpc.beq_0(:) = mpc.beq_0(:) - mpc.Bd(:,:,1)*mpc.d(:,1);

    for k = 1:mpc.N-1
        mpc.beq_k(:,k) = - mpc.Bd(:,:,k+1)*mpc.d(:,k+1);
    end
end

%% inequalities
% k = 0

if mpc.has_du_cnstr
    if mpc.du_cnstr.min_limit
        row = mpc.du_cnstr.min_ineqRow_0;
        %Ineq [-I I]*[u0 g0]'=-du_min-u_prev
        mpc.bi_0(row) = -mpc.du_cnstr.min(:,1)-u_prev;
    end
    if mpc.du_cnstr.max_limit
        row = mpc.du_cnstr.max_ineqRow_0;
        %Ineq [I I]*[u0 g0]'=du_max+u_prev
        mpc.bi_0(row) = mpc.du_cnstr.max(:,1)+u_prev;
    end
end

if mpc.has_y_cnstr

    if mpc.y_cnstr.min_limit && mpc.y_use_k0
        row = mpc.y_cnstr.min_ineqRow_0;
        %Ineq [-D I -I]*[u g v]' = -y_min +C*s+Dd*d
        mpc.bi_0(row) = -mpc.y_cnstr.min_0;
        if mpc.y_use_s, mpc.bi_0(row) = mpc.bi_0(row) + mpc.C_0*s_prev; end
        if mpc.y_use_d, mpc.bi_0(row) = mpc.bi_0(row) + mpc.Dd_0*mpc.d(:,1); end
    end

    if mpc.y_cnstr.max_limit && mpc.y_use_k0
        row = mpc.y_cnstr.max_ineqRow_0;
        %Ineq [D I -I]*[u g v]' = y_max -C*s-Dd*d
        mpc.bi_0(row) = mpc.y_cnstr.max_0;
        if mpc.y_use_s, mpc.bi_0(row) = mpc.bi_0(row) - mpc.C_0*s_prev; end
        if mpc.y_use_d, mpc.bi_0(row) = mpc.bi_0(row) - mpc.Dd_0*mpc.d(:,1); end
    end
end

if mpc.has_h_cnstr

    if mpc.h_cnstr.min_limit && mpc.h_cnstr.use_k0
        row = mpc.h_cnstr.min_ineqRow_0;
        %Ineq [-D_0 I -I]*[u g v]' = -h_min_0 +C_0*s+Dsu_0*su+Dd_0*d
        mpc.bi_0(row) = -mpc.h_cnstr.min_0;
        if mpc.h_cnstr.use_s, mpc.bi_0(row) = mpc.bi_0(row) + mpc.Ch_0*s_prev; end
        if mpc.h_cnstr.use_su, mpc.bi_0(row) = mpc.bi_0(row) + mpc.Dsuh_0*u_prev; end
        if mpc.h_cnstr.use_d, mpc.bi_0(row) = mpc.bi_0(row) + mpc.Ddh_0*mpc.dh(:,1); end
    end

    if mpc.h_cnstr.max_limit && mpc.h_cnstr.use_k0
        row = mpc.h_cnstr.max_ineqRow_0;
        %Ineq [D_0 I -I]*[u g v]' = h_max_0 -C_0*s-Dsu_0*su-Dd_0*d
        mpc.bi_0(row) = mpc.h_cnstr.max_0;
        if mpc.h_cnstr.use_s, mpc.bi_0(row) = mpc.bi_0(row) - mpc.Ch_0*s_prev; end
        if mpc.h_cnstr.use_su, mpc.bi_0(row) = mpc.bi_0(row) - mpc.Dsuh_0*u_prev; end
        if mpc.h_cnstr.use_d, mpc.bi_0(row) = mpc.bi_0(row) - mpc.Ddh_0*mpc.dh(:,1); end
    end
end

for k = 1:mpc.N-1

    if mpc.has_h_cnstr

        if mpc.h_cnstr.min_limit && mpc.h_cnstr.use_d
            row = mpc.h_cnstr.min_ineqRow_k;
            %Ineq [-C -Dsu -D I -I]*[s su u g v]' = -h_min+Dd*d
            mpc.bi_k(row,k) = -mpc.h_cnstr.min(:,k) + mpc.Ddh(:,:,k)*mpc.dh(:,k+1);
        end

        if mpc.h_cnstr.max_limit && mpc.h_cnstr.use_d
            row = mpc.h_cnstr.max_ineqRow_k;
            %Ineq [C Dsu D I -I]*[s su u g v]' = h_max-Dd*d
            mpc.bi_k(row,k) = mpc.h_cnstr.max(:,k) - mpc.Ddh(:,:,k)*mpc.dh(:,k+1);
        end
    end

    if mpc.has_y_cnstr

        if mpc.y_cnstr.min_limit && mpc.y_use_d
            row = mpc.y_cnstr.min_ineqRow_k;
            %Ineq [-C -D I -I]*[s u g v]' = -y_min+Dd*d
            mpc.bi_k(row,k) = -mpc.y_cnstr.min(:,k) + mpc.Dd(:,:,k)*mpc.d(:,k+1);
        end

        if mpc.y_cnstr.max_limit && mpc.y_use_d
            row = mpc.y_cnstr.max_ineqRow_k;
            %Ineq [C D I -I]*[s u g v]' = y_max-Dd*d
            mpc.bi_k(row,k) = mpc.y_cnstr.max(:,k) - mpc.Dd(:,:,k)*mpc.d(:,k+1);
        end
    end
end

end