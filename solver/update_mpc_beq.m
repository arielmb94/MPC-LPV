function mpc = update_mpc_beq(mpc,s_prev,u_prev)

% k = 0
index = mpc.dyn_k(:,1);
mpc.beq(index) = -mpc.A*s_prev;
if mpc.dyn_use_d
    mpc.beq(index) = mpc.beq(index) - mpc.Bd*mpc.d(:,1);
end

if mpc.has_du_cnstr
    if mpc.du_cnstr.min_limit
        index = mpc.du_cnstr.min_eq_index_k(:,1);
        mpc.beq(index) = -mpc.du_cnstr.min-u_prev;
    end
    if mpc.du_cnstr.max_limit
        index = mpc.du_cnstr.max_eq_index_k(:,1);
        mpc.beq(index) = mpc.du_cnstr.max+u_prev;
    end
end

if mpc.has_h_cnstr

    if mpc.h_cnstr.min_limit && mpc.h_cnstr.use_u
        index = mpc.h_cnstr.min_eq_index_k(:,1);
        mpc.beq(index) = -mpc.h_cnstr.min;
        if mpc.h_cnstr.use_s
            mpc.beq(index) = mpc.beq(index) + mpc.Ch*s_prev;
        end
        if mpc.h_cnstr.use_su
            mpc.beq(index) = mpc.beq(index) + mpc.Dduh*u_prev;
        end
        if mpc.h_cnstr.use_d
            mpc.beq(index) = mpc.beq(index) + mpc.Ddh*mpc.dh(:,1);
        end
    end

    if mpc.h_cnstr.max_limit && mpc.h_cnstr.use_u
        index = mpc.h_cnstr.max_eq_index_k(:,1);
        mpc.beq(index) = mpc.h_cnstr.max;
        if mpc.h_cnstr.use_s
            mpc.beq(index) = mpc.beq(index) - mpc.Ch*s_prev;
        end
        if mpc.h_cnstr.use_su
            mpc.beq(index) = mpc.beq(index) - mpc.Dduh*u_prev;
        end
        if mpc.h_cnstr.use_d
            mpc.beq(index) = mpc.beq(index) - mpc.Ddh*mpc.dh(:,1);
        end
    end
end

for k = 2:mpc.N

    if mpc.dyn_use_d
        index = mpc.dyn_k(:,k);
        mpc.beq(index) = - mpc.Bd*mpc.d(:,k);
    end

    if mpc.has_h_cnstr

        if mpc.h_cnstr.min_limit && mpc.h_cnstr.use_d
            index = mpc.h_cnstr.min_eq_index_k(:,k);
            mpc.beq(index) = -mpc.h_cnstr.min + mpc.Ddh*mpc.dh(:,k);
        end

        if mpc.h_cnstr.max_limit && mpc.h_cnstr.use_d
            index = mpc.h_cnstr.max_eq_index_k(:,k);
            mpc.beq(index) = mpc.h_cnstr.max - mpc.Ddh*mpc.dh(:,k);
        end
    end

    if mpc.has_y_cnstr

        if mpc.y_cnstr.min_limit && mpc.y_cnstr.use_d
            index = mpc.y_cnstr.min_eq_index_k(:,k);
            mpc.beq(index) = -mpc.y_cnstr.min + mpc.Dd*mpc.d(:,k);
        end

        if mpc.y_cnstr.max_limit && mpc.y_cnstr.use_d
            index = mpc.y_cnstr.max_eq_index_k(:,k);
            mpc.beq(index) = mpc.y_cnstr.max - mpc.Dd*mpc.d(:,k);
        end
    end

end

end