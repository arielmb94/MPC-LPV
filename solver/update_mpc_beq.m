function mpc = update_mpc_beq(mpc,s_prev,u_prev)

% k = 0
mpc.beq(1:mpc.nx,1) = -mpc.A*s_prev;

if mpc.dyn_use_d
    mpc.beq(1:mpc.nx,1) = mpc.beq(1:mpc.nx,1) - mpc.Bd*mpc.d(:,1);

    for k = 2:mpc.N
        mpc.beq(1:mpc.nx,k) = - mpc.Bd*mpc.d(:,k);
    end
end

%%
% k = 0

if mpc.has_du_cnstr
    if mpc.du_cnstr.min_limit
        row = mpc.du_cnstr.min_row_0;
        %Ineq [-I I]*[u0 g0]'=-du_min-u_prev
        mpc.bi_0(row) = -mpc.du_cnstr.min-u_prev;
    end
    if mpc.du_cnstr.max_limit
        row = mpc.du_cnstr.max_row_0;
        %Ineq [I I]*[u0 g0]'=du_max+u_prev
        mpc.bi_0(row) = mpc.du_cnstr.max+u_prev;
    end
end

if mpc.has_h_cnstr

    if mpc.h_cnstr.min_limit && mpc.h_cnstr.use_u
        row = mpc.h_cnstr.min_row_0;
        %Ineq [-D I -I]*[u g v]' = -h_min +Cs+Dsu*su+Dd*d
        mpc.bi_0(row) = -mpc.h_cnstr.min;
        if mpc.h_cnstr.use_s
            mpc.bi_0(row) = mpc.bi_0(row) + mpc.Ch*s_prev;
        end
        if mpc.h_cnstr.use_su
            mpc.bi_0(row) = mpc.bi_0(row) + mpc.Dduh*u_prev;
        end
        if mpc.h_cnstr.use_d
            mpc.bi_0(row) = mpc.bi_0(row) + mpc.Ddh*mpc.dh(:,1);
        end
    end

    if mpc.h_cnstr.max_limit && mpc.h_cnstr.use_u
        row = mpc.h_cnstr.max_row_0;
        %Ineq [D I -I]*[u g v]' = h_max -Cs-Dsu*su-Dd*d
        mpc.bi_0(row) = mpc.h_cnstr.max;
        if mpc.h_cnstr.use_s
            mpc.bi_0(row) = mpc.bi_0(row) - mpc.Ch*s_prev;
        end
        if mpc.h_cnstr.use_su
            mpc.bi_0(row) = mpc.bi_0(row) - mpc.Dduh*u_prev;
        end
        if mpc.h_cnstr.use_d
            mpc.bi_0(row) = mpc.bi_0(row) - mpc.Ddh*mpc.dh(:,1);
        end
    end
end

for k = 2:mpc.N

    if mpc.has_h_cnstr

        if mpc.h_cnstr.min_limit && mpc.h_cnstr.use_d
            row = mpc.h_cnstr.min_row_k;
            %Ineq [-C -Dsu -D I -I]*[s su u g v]' = -h_min+Dd*d
            mpc.bi_k(row,k-1) = -mpc.h_cnstr.min + mpc.Ddh*mpc.dh(:,k);
        end

        if mpc.h_cnstr.max_limit && mpc.h_cnstr.use_d
            row = mpc.h_cnstr.max_row_k;
            %Ineq [C Dsu D I -I]*[s su u g v]' = h_max-Dd*d
            mpc.bi_k(row,k-1) = mpc.h_cnstr.max - mpc.Ddh*mpc.dh(:,k);
        end
    end

    if mpc.has_y_cnstr

        if mpc.y_cnstr.min_limit && mpc.y_cnstr.use_d
            row = mpc.y_cnstr.min_row_k;
            %Ineq [-C -D I -I]*[s u g v]' = -y_min+Dd*d
            mpc.bi_k(row,k-1) = -mpc.y_cnstr.min + mpc.Dd*mpc.d(:,k);
        end

        if mpc.y_cnstr.max_limit && mpc.y_cnstr.use_d
            row = mpc.y_cnstr.max_row_k;
            %Ineq [C D I -I]*[s u g v]' = y_max-Dd*d
            mpc.bi_k(row,k-1) = mpc.y_cnstr.max - mpc.Dd*mpc.d(:,k);
        end
    end

end

end