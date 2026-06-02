function mpc = genEqualities(mpc)

du = mpc.has_du;
nx = mpc.nx;
nu = mpc.nu;
ny = mpc.ny;
nh = mpc.nh;
N = mpc.N;

r_len = (mpc.Nx) + (mpc.Nu-1)*du + sum(mpc.ng_k);
c_len = mpc.Nx+mpc.Nu+(mpc.Nu-mpc.nu)*du+sum(mpc.ng_k)+sum(mpc.nv_k);

Aeq = zeros(r_len,c_len);
beq = zeros(r_len,1);

% k = 0
%Dynamics [B -I]*[u0 s1]'=-A*s0-Dd*d0
row = mpc.dyn_k(:,1);
u_col = mpc.u_index_k(:,1);
s_next_col = mpc.s_index_k(:,2);

Aeq = appendDynamics(Aeq, row,...
                        [], u_col, s_next_col,...
                        [], mpc.B, -eye(nx));

if mpc.has_du
    %[I 0 0 -I]*[u0 g0 s1 su1]'
    row = mpc.su_dyn_k(:,1);
    u_col = mpc.u_index_k(:,1);
    s_next_col = mpc.su_index_k(:,2);

    Aeq = appendDynamics(Aeq, row,...
                        [], u_col, s_next_col,...
                        [], eye(nu), -eye(nu));
end

if mpc.has_u_cnstr
    if mpc.u_cnstr.min_limit
        %Ineq [-I I]*[u0 g0]'=-u_min
        row = mpc.u_cnstr.min_eq_index_k(:,1);
        u_col = mpc.u_index_k(:,1);
        g_col = mpc.u_cnstr.g_min_index_k(:,1);
        b_val = -mpc.u_cnstr.min;

        [Aeq, beq] = appendGeneralizedConstraint(Aeq, beq, row,...
            [], u_col, [], g_col, [],...
            [], -eye(nu), [], eye(nu), [], b_val);

    end
    if mpc.u_cnstr.max_limit
        %Ineq [I I]*[u0 g0]'=u_max
        row = mpc.u_cnstr.max_eq_index_k(:,1);
        u_col = mpc.u_index_k(:,1);
        g_col = mpc.u_cnstr.g_max_index_k(:,1);
        b_val = mpc.u_cnstr.max;

        [Aeq, beq] = appendGeneralizedConstraint(Aeq, beq, row,...
            [], u_col, [], g_col, [],...
            [], eye(nu), [], eye(nu), [], b_val);
    end
end

if mpc.has_du_cnstr
    if mpc.du_cnstr.min_limit
        %Ineq [-I I]*[u0 g0]'=-du_min-u_prev
        row = mpc.du_cnstr.min_eq_index_k(:,1);
        u_col = mpc.u_index_k(:,1);
        g_col = mpc.du_cnstr.g_min_index_k(:,1);
        b_val = -mpc.du_cnstr.min;

        [Aeq, beq] = appendGeneralizedConstraint(Aeq, beq, row,...
            [], u_col, [], g_col, [],...
            [], -eye(nu), [], eye(nu), [], b_val);
    end
    if mpc.du_cnstr.max_limit
        %Ineq [I I]*[u0 g0]'=du_max+u_prev
        row = mpc.du_cnstr.max_eq_index_k(:,1);
        u_col = mpc.u_index_k(:,1);
        g_col = mpc.du_cnstr.g_max_index_k(:,1);
        b_val = mpc.du_cnstr.max;

        [Aeq, beq] = appendGeneralizedConstraint(Aeq, beq, row,...
            [], u_col, [], g_col, [],...
            [], eye(nu), [], eye(nu), [], b_val);
    end
end

if mpc.has_h_cnstr && any(mpc.Dh)
    if mpc.h_cnstr.min_limit
        %Ineq [-D I -I]*[u g v]' = -h_min +Cs+Ddu*su+Dd*d
        row = mpc.h_cnstr.min_eq_index_k(:,1);
        u_col = mpc.u_index_k(:,1);
        g_col = mpc.h_cnstr.g_min_index_k(:,1);
        v_col = mpc.h_cnstr.v_min_index_k(:,1);
        b_val = -mpc.h_cnstr.min;

        [Aeq, beq] = appendGeneralizedConstraint(Aeq, beq, row,...
            [], u_col, [], g_col, v_col,...
            [],-mpc.Dh,[],eye(nh),-eye(nh), b_val);
    end
    if mpc.h_cnstr.max_limit
        %Ineq [D I -I]*[u g v]' = h_max -Cs-Ddu*su-Dd*d
        row = mpc.h_cnstr.max_eq_index_k(:,1);
        u_col = mpc.u_index_k(:,1);
        g_col = mpc.h_cnstr.g_max_index_k(:,1);
        v_col = mpc.h_cnstr.v_max_index_k(:,1);
        b_val = mpc.h_cnstr.max;

        [Aeq, beq] = appendGeneralizedConstraint(Aeq, beq, row,...
            [], u_col, [], g_col, v_col,...
            [], mpc.Dh,[],eye(nh),-eye(nh), b_val);
    end
end


for k = 2:N
    %Dynamics [A B -I]*[s1 u1 s2]'
    row = mpc.dyn_k(:,k);
    s_col = mpc.s_index_k(:,k);
    u_col = mpc.u_index_k(:,k);
    s_next_col = mpc.s_index_k(:,k+1);

    Aeq = appendDynamics(Aeq, row,...
        s_col, u_col, s_next_col,...
        mpc.A, mpc.B, -eye(nx));

    %[I -I]*[uk su1]'
    if mpc.has_du && all(mpc.su_index_k(:,k+1))
        row = mpc.su_dyn_k(:,k);
        u_col = mpc.u_index_k(:,k);
        s_next_col = mpc.su_index_k(:,k+1);

        Aeq = appendDynamics(Aeq, row,...
            [], u_col, s_next_col,...
            [], eye(nu), -eye(nu));
    end

    if mpc.has_s_cnstr
        if mpc.s_cnstr.min_limit
            %Ineq [-I I -I]*[sk gk vk]'=-s_min
            row = mpc.s_cnstr.min_eq_index_k(:,k);
            s_col = mpc.s_index_k(:,k);
            g_col = mpc.s_cnstr.g_min_index_k(:,k);
            v_col = mpc.s_cnstr.v_min_index_k(:,k);
            b_val = -mpc.s_cnstr.min;

            [Aeq, beq] = appendGeneralizedConstraint(Aeq, beq, row,...
                s_col, [], [], g_col, v_col,...
                -eye(nx), [], [], eye(nx), -eye(nx), b_val);
        end
        if mpc.s_cnstr.max_limit
            %Ineq [I I -I]*[sk gk vk]'=s_max
            row = mpc.s_cnstr.max_eq_index_k(:,k);
            s_col = mpc.s_index_k(:,k);
            g_col = mpc.s_cnstr.g_max_index_k(:,k);
            v_col = mpc.s_cnstr.v_max_index_k(:,k);
            b_val = mpc.s_cnstr.max;

            [Aeq, beq] = appendGeneralizedConstraint(Aeq, beq, row,...
                s_col, [], [], g_col, v_col,...
                eye(nx), [], [], eye(nx), -eye(nx), b_val);
        end
    end

    if mpc.has_u_cnstr
        if mpc.u_cnstr.min_limit
            %Ineq [-I I]*[u0 g0]'=-u_min
            row = mpc.u_cnstr.min_eq_index_k(:,k);
            u_col = mpc.u_index_k(:,k);
            g_col = mpc.u_cnstr.g_min_index_k(:,k);
            b_val = -mpc.u_cnstr.min;

            [Aeq, beq] = appendGeneralizedConstraint(Aeq, beq, row,...
                [], u_col, [], g_col, [],...
                [], -eye(nu), [], eye(nu), [], b_val);
        end
        if mpc.u_cnstr.max_limit
            %Ineq [I I]*[u0 g0]'=u_max
            row = mpc.u_cnstr.max_eq_index_k(:,k);
            u_col = mpc.u_index_k(:,k);
            g_col = mpc.u_cnstr.g_max_index_k(:,k);
            b_val = mpc.u_cnstr.max;

            [Aeq, beq] = appendGeneralizedConstraint(Aeq, beq, row,...
                [], u_col, [], g_col, [],...
                [], eye(nu), [], eye(nu), [], b_val);
        end
    end

    if mpc.has_du_cnstr
        if mpc.du_cnstr.min_limit
            %Ineq [I -I I]*[su u g]' = -du_min
            row = mpc.du_cnstr.min_eq_index_k(:,k);
            u_col = mpc.u_index_k(:,k);
            su_col = mpc.su_index_k(:,k);
            g_col = mpc.du_cnstr.g_min_index_k(:,k);
            b_val = -mpc.du_cnstr.min;

            [Aeq, beq] = appendGeneralizedConstraint(Aeq, beq, row,...
                [], u_col, su_col, g_col, [],...
                [],-eye(nu),eye(nu),eye(nu),[], b_val);
        end
        if mpc.du_cnstr.max_limit
            %Ineq [-I I I]*[su u g]' = du_max
            row = mpc.du_cnstr.max_eq_index_k(:,k);
            u_col = mpc.u_index_k(:,k);
            su_col = mpc.su_index_k(:,k);
            g_col = mpc.du_cnstr.g_max_index_k(:,k);
            b_val = mpc.du_cnstr.max;

            [Aeq, beq] = appendGeneralizedConstraint(Aeq, beq, row,...
                [], u_col, su_col, g_col, [],...
                [],eye(nu),-eye(nu),eye(nu),[], b_val);
        end
    end

    if mpc.has_y_cnstr
        if mpc.y_cnstr.min_limit
            %Ineq [-C -D I -I]*[s u g v]' = -y_min+Dd*d
            row = mpc.y_cnstr.min_eq_index_k(:,k);
            s_col = mpc.s_index_k(:,k);
            u_col = mpc.u_index_k(:,k);
            g_col = mpc.y_cnstr.g_min_index_k(:,k);
            v_col = mpc.y_cnstr.v_min_index_k(:,k);
            b_val = -mpc.y_cnstr.min;

            [Aeq, beq] = appendGeneralizedConstraint(Aeq, beq, row,...
                s_col, u_col, [], g_col, v_col,...
                -mpc.C,-mpc.D,[],eye(ny),-eye(ny), b_val);
        end
        if mpc.y_cnstr.max_limit
            %Ineq [C D I -I]*[s u g v]' = y_max-Dd*d
            row = mpc.y_cnstr.max_eq_index_k(:,k);
            s_col = mpc.s_index_k(:,k);
            u_col = mpc.u_index_k(:,k);
            g_col = mpc.y_cnstr.g_max_index_k(:,k);
            v_col = mpc.y_cnstr.v_max_index_k(:,k);
            b_val = mpc.y_cnstr.max;

            [Aeq, beq] = appendGeneralizedConstraint(Aeq, beq, row,...
                s_col, u_col, [], g_col, v_col,...
                mpc.C,mpc.D,[],eye(ny),-eye(ny), b_val);
        end
    end

    if mpc.has_h_cnstr
        if mpc.h_cnstr.min_limit
            %Ineq [-C -Ddu -D I -I]*[s su u g v]' = -h_min+Dd*d
            row = mpc.h_cnstr.min_eq_index_k(:,k);
            s_col = mpc.s_index_k(:,k);
            su_col = mpc.su_index_k(:,k);
            u_col = mpc.u_index_k(:,k);
            g_col = mpc.h_cnstr.g_min_index_k(:,k);
            v_col = mpc.h_cnstr.v_min_index_k(:,k);
            b_val = -mpc.h_cnstr.min;

            [Aeq, beq] = appendGeneralizedConstraint(Aeq, beq, row,...
                s_col, u_col, su_col, g_col, v_col,...
                -mpc.Ch,-mpc.Dh,-mpc.Dduh,eye(nh),-eye(nh), b_val);
        end
        if mpc.h_cnstr.max_limit
            %Ineq [C Ddu D I -I]*[s su u g v]' = y_max-Dd*d
            row = mpc.h_cnstr.max_eq_index_k(:,k);
            s_col = mpc.s_index_k(:,k);
            su_col = mpc.su_index_k(:,k);
            u_col = mpc.u_index_k(:,k);
            g_col = mpc.h_cnstr.g_max_index_k(:,k);
            v_col = mpc.h_cnstr.v_max_index_k(:,k);
            b_val = mpc.h_cnstr.max;

            [Aeq, beq] = appendGeneralizedConstraint(Aeq, beq, row,...
                s_col, u_col, su_col, g_col, v_col,...
                mpc.Ch,mpc.Dh,mpc.Dduh,eye(nh),-eye(nh), b_val);
        end
    end

end

if mpc.has_s_cnstr
    if mpc.s_cnstr.min_limit
        %Ineq [-I I -I]*[sk gk vk]'=-s_min
        row = mpc.s_cnstr.min_eq_index_k(:,N+1);
        s_col = mpc.s_index_k(:,N+1);
        g_col = mpc.s_cnstr.g_min_index_k(:,N+1);
        v_col = mpc.s_cnstr.v_min_index_k(:,N+1);
        b_val = -mpc.s_cnstr.min;

        [Aeq, beq] = appendGeneralizedConstraint(Aeq, beq, row,...
            s_col, [], [], g_col, v_col,...
            -eye(nx), [], [], eye(nx), -eye(nx), b_val);
    end
    if mpc.s_cnstr.max_limit
        %Ineq [I I -I]*[sk gk vk]'=s_max
        row = mpc.s_cnstr.max_eq_index_k(:,N+1);
        s_col = mpc.s_index_k(:,N+1);
        g_col = mpc.s_cnstr.g_max_index_k(:,N+1);
        v_col = mpc.s_cnstr.v_max_index_k(:,N+1);
        b_val = mpc.s_cnstr.max;

        [Aeq, beq] = appendGeneralizedConstraint(Aeq, beq, row,...
            s_col, [], [], g_col, v_col,...
            eye(nx), [], [], eye(nx), -eye(nx), b_val);
    end
end

mpc.Aeq = Aeq;
mpc.beq = beq;
    
end



