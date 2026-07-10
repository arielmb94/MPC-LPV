function mpc = genEqualities(mpc)

du = mpc.has_du;
nx = mpc.nx;
nu = mpc.nu;
nse = nx+nu*du;
ny = mpc.ny;
nh = mpc.nh;
N = mpc.N;

mpc.nse = nse;

r_len = (mpc.Nx) + (mpc.Nu)*du + sum(mpc.ng_k);
c_len = mpc.Nx+mpc.Nu+(mpc.Nu)*du+sum(mpc.ng_k)+sum(mpc.nv_k);

beq = zeros(nse,mpc.N);
bi_0 = zeros(mpc.ng_k(1),1);
bi_k = zeros(mpc.ng_k(2),mpc.N-1);

mpc.rp_0 = zeros(nse,1);
mpc.rp_k = zeros(nse,mpc.N-1);
mpc.beq_0 = zeros(nx,1);
mpc.beq_k = zeros(nx,mpc.N-1);
mpc.ri_0 = zeros(mpc.ng_k(1),1);
mpc.ri_k = zeros(mpc.ng_k(2),mpc.N-1);
mpc.ri_ter = zeros(mpc.ng_k(mpc.N+1),1);

mpc.ri_hat_0 = zeros(mpc.ng_k(1),1);
mpc.ri_hat_k = zeros(mpc.ng_k(2),mpc.N-1);
mpc.ri_hat_ter = zeros(mpc.ng_k(mpc.N+1),1);

mpc.S_0 = zeros(mpc.ng_k(1),1);
mpc.S_k = zeros(mpc.ng_k(2),mpc.N-1);
mpc.S_ter = zeros(mpc.ng_k(mpc.N+1),1);

mpc.iS_0 = zeros(mpc.ng_k(1),1);
mpc.iS_k = zeros(mpc.ng_k(2),mpc.N-1);
mpc.iS_ter = zeros(mpc.ng_k(mpc.N+1),1);

vi_0 = [];

%k = 0
Ai_0 = zeros(mpc.ng_k(1),mpc.nu);
start_index = 1;

u_col = 1:mpc.nu;

if mpc.has_u_cnstr
    if mpc.u_cnstr.min_limit
        %Ineq [-I I]*[u0 g0]'=-u_min
        row = start_index:start_index+mpc.nu-1;

        gi_index_k = mpc.g_index_0(row);
        mpc.u_cnstr.g_min_index_k = gi_index_k;
        mpc.u_cnstr.min_row_0 = row;

        b_val = -mpc.u_cnstr.min;

        [Ai_0, bi_0] = appendGeneralizedConstraint(Ai_0, bi_0, row,...  
                                                 [], [], u_col, ...
                                                 [], [], -eye(nu), b_val);

        start_index = start_index + nu;
    end
    if mpc.u_cnstr.max_limit
        %Ineq [I I]*[u0 g0]'=u_max
        row = start_index:start_index+mpc.nu-1;

        gi_index_k = mpc.g_index_0(row);
        mpc.u_cnstr.g_max_index_k = gi_index_k;
        mpc.u_cnstr.max_row_0 = row;

        b_val = mpc.u_cnstr.max;

        [Ai_0, bi_0] = appendGeneralizedConstraint(Ai_0, bi_0, row,...  
                                                 [], [], u_col, ...
                                                 [], [], eye(nu), b_val);

        start_index = start_index + nu;
    end
end

if mpc.has_du_cnstr
    if mpc.du_cnstr.min_limit
        %Ineq [-I I]*[u0 g0]'=-du_min-u_prev
        row = start_index:start_index+mpc.nu-1;

        gi_index_k = mpc.g_index_0(row);
        mpc.du_cnstr.g_min_index_k = gi_index_k;
        mpc.du_cnstr.min_row_0 = row;

        b_val = -mpc.du_cnstr.min;

        [Ai_0, bi_0] = appendGeneralizedConstraint(Ai_0, bi_0, row,...  
                                                 [], [], u_col, ...
                                                 [], [], -eye(nu), b_val);

        start_index = start_index + nu;
    end
    if mpc.du_cnstr.max_limit
        %Ineq [I I]*[u0 g0]'=du_max+u_prev
        row = start_index:start_index+mpc.nu-1;

        gi_index_k = mpc.g_index_0(row);
        mpc.du_cnstr.g_max_index_k = gi_index_k;
        mpc.du_cnstr.max_row_0 = row;

        b_val = mpc.du_cnstr.max;

        [Ai_0, bi_0] = appendGeneralizedConstraint(Ai_0, bi_0, row,...  
                                                 [], [], u_col, ...
                                                 [], [], eye(nu), b_val);

        start_index = start_index + nu;
    end
end

start_index_v = 1;
if mpc.has_h_cnstr && any(mpc.Dh)
    if mpc.h_cnstr.min_limit
        %Ineq [-D I -I]*[u g v]' = -h_min +Cs+Dsu*su+Dd*d
        row = start_index:start_index+mpc.nh-1;
        row_v = start_index_v:start_index_v+mpc.nh-1;

        mpc.h_cnstr.min_row_0 = row;
        mpc.h_cnstr.min_row_v_0 = row_v;

        gi_index_k = mpc.g_index_0(row);
        mpc.h_cnstr.g_min_index_k = gi_index_k;

        vi_index_k = mpc.v_index_0(row_v);
        mpc.h_cnstr.v_min_index_k = vi_index_k;
        vi_0 = [vi_0; row];

        D_mat = -mpc.Dh;
        b_val = -mpc.h_cnstr.min;

        [Ai_0, bi_0] = appendGeneralizedConstraint(Ai_0, bi_0, row,...  
                                                 [], [], u_col, ...
                                                 [], [], D_mat, b_val);
        start_index = start_index + nh;
        start_index_v = start_index_v + nh;
    end
    if mpc.h_cnstr.max_limit
        %Ineq [D I -I]*[u g v]' = h_max -Cs-Dsu*su-Dd*d
        row = start_index:start_index+mpc.nh-1;
        row_v = start_index_v:start_index_v+mpc.nh-1;

        mpc.h_cnstr.max_row_0 = row;
        mpc.h_cnstr.max_row_v_0 = row_v;

        gi_index_k = mpc.g_index_0(row);
        mpc.h_cnstr.g_max_index_k = gi_index_k;

        vi_index_k = mpc.v_index_0(row_v);
        mpc.h_cnstr.v_max_index_k = vi_index_k;
        vi_0 = [vi_0; row];

        D_mat = mpc.Dh;
        b_val = mpc.h_cnstr.max;

        [Ai_0, bi_0] = appendGeneralizedConstraint(Ai_0, bi_0, row,...  
                                                 [], [], u_col, ...
                                                 [], [], D_mat, b_val);
        start_index = start_index + nh;
        start_index_v = start_index_v + nh;
    end
end
mpc.vi_0 = vi_0;
mpc.Ai_0 = Ai_0;  
mpc.bi_0 = bi_0;

s_col = 1:nx;
su_col = nx+1:nse;
u_col = nse+1:nse+nu;

mpc.inq_s_col = s_col;
mpc.inq_su_col = su_col;
mpc.inq_u_col = u_col;

Ai_k = [];
for k = 1:N-1

    Ai = zeros(mpc.ng_k(2),nx+nu*du+nu);
    bi = zeros(mpc.ng_k(2),1);

    vi_k = [];
    start_index = 1;
    start_index_v = 1;

    if mpc.has_s_cnstr
        if mpc.s_cnstr.min_limit
            %Ineq [-I I -I]*[sk gk vk]'=-s_min
            row = start_index:start_index+mpc.nx-1;
            row_v = start_index_v:start_index_v+mpc.nx-1;

            mpc.s_cnstr.min_row_k = row;
            mpc.s_cnstr.min_row_v_k = row_v;
           
            gi_index_k = mpc.g_index_k(row,k);
            mpc.s_cnstr.g_min_index_k = [mpc.s_cnstr.g_min_index_k gi_index_k];

            vi_index_k = mpc.v_index_k(row_v,k);
            mpc.s_cnstr.v_min_index_k = [mpc.s_cnstr.v_min_index_k vi_index_k];
            vi_k = [vi_k; row'];

            b_val = -mpc.s_cnstr.min;

            [Ai, bi] = appendGeneralizedConstraint(Ai, bi, row,...  
                                                   s_col, [], [], ...
                                                   -eye(nx), [], [], b_val);
            start_index = start_index + nx;
            start_index_v = start_index_v + nx;
        end
        if mpc.s_cnstr.max_limit
            %Ineq [I I -I]*[sk gk vk]'=s_max
            row = start_index:start_index+mpc.nx-1;
            row_v = start_index_v:start_index_v+mpc.nx-1;

            mpc.s_cnstr.max_row_k = row;
            mpc.s_cnstr.max_row_v_k = row_v;

            gi_index_k = mpc.g_index_k(row,k);
            mpc.s_cnstr.g_max_index_k = [mpc.s_cnstr.g_max_index_k gi_index_k];

            vi_index_k = mpc.v_index_k(row_v,k);
            mpc.s_cnstr.v_max_index_k = [mpc.s_cnstr.v_max_index_k vi_index_k];
            vi_k = [vi_k; row'];

            b_val = mpc.s_cnstr.max;

            [Ai, bi] = appendGeneralizedConstraint(Ai, bi, row,...  
                                                   s_col, [], [], ...
                                                   eye(nx), [], [], b_val);
            start_index = start_index + nx;
            start_index_v = start_index_v + nx;
        end
    end

    if mpc.has_u_cnstr
        if mpc.u_cnstr.min_limit
            %Ineq [-I I]*[u0 g0]'=-u_min
            row = start_index:start_index+mpc.nu-1;

            mpc.u_cnstr.min_row_k = row;

            gi_index_k = mpc.g_index_k(row,k);
            mpc.u_cnstr.g_min_index_k = [mpc.u_cnstr.g_min_index_k gi_index_k];

            b_val = -mpc.u_cnstr.min;

            [Ai, bi] = appendGeneralizedConstraint(Ai, bi, row,...
                [], [], u_col, ...
                [], [], -eye(nu), b_val);

            start_index = start_index + nu;
        end
        if mpc.u_cnstr.max_limit
            %Ineq [I I]*[u0 g0]'=u_max
            row = start_index:start_index+mpc.nu-1;

            mpc.u_cnstr.max_row_k = row;

            gi_index_k = mpc.g_index_k(row,k);
            mpc.u_cnstr.g_max_index_k = [mpc.u_cnstr.g_max_index_k gi_index_k];

            b_val = mpc.u_cnstr.max;

            [Ai, bi] = appendGeneralizedConstraint(Ai, bi, row,...
                [], [], u_col, ...
                [], [], eye(nu), b_val);

            start_index = start_index + nu;
        end
    end

    if mpc.has_du_cnstr
        if mpc.du_cnstr.min_limit
            %Ineq [I -I I]*[su u g]' = -du_min
            row = start_index:start_index+mpc.nu-1;

            mpc.du_cnstr.min_row_k = row;

            gi_index_k = mpc.g_index_k(row,k);
            mpc.du_cnstr.g_min_index_k = [mpc.du_cnstr.g_min_index_k gi_index_k];

            b_val = -mpc.du_cnstr.min;

            [Ai, bi] = appendGeneralizedConstraint(Ai, bi, row,...
                [], su_col, u_col, ...
                [], eye(nu), -eye(nu), b_val);

            start_index = start_index + nu;
        end
        if mpc.du_cnstr.max_limit
            %Ineq [-I I I]*[su u g]' = du_max
            row = start_index:start_index+mpc.nu-1;

            mpc.du_cnstr.max_row_k = row;

            gi_index_k = mpc.g_index_k(row,k);
            mpc.du_cnstr.g_max_index_k = [mpc.du_cnstr.g_max_index_k gi_index_k];

            b_val = mpc.du_cnstr.max;

            [Ai, bi] = appendGeneralizedConstraint(Ai, bi, row,...
                [], su_col, u_col, ...
                [], -eye(nu), eye(nu), b_val);

            start_index = start_index + nu;
        end
    end

    if mpc.has_y_cnstr
        if mpc.y_cnstr.min_limit
            %Ineq [-C -D I -I]*[s u g v]' = -y_min+Dd*d
            row = start_index:start_index+ny-1;
            row_v = start_index_v:start_index_v+mpc.ny-1;

            mpc.y_cnstr.min_row_k = row;
            mpc.y_cnstr.min_row_v_k = row_v;

            gi_index_k = mpc.g_index_k(row,k);
            mpc.y_cnstr.g_min_index_k = [mpc.y_cnstr.g_min_index_k gi_index_k];
            mpc.y_cnstr.ineq_min_row = row;

            vi_index_k = mpc.v_index_k(row_v,k);
            mpc.y_cnstr.v_min_index_k = [mpc.y_cnstr.v_min_index_k vi_index_k];
            vi_k = [vi_k; row'];

            C_mat = -mpc.C;
            D_mat = -mpc.D;
            b_val = -mpc.y_cnstr.min;

            [Ai, bi] = appendGeneralizedConstraint(Ai, bi, row,...
                                                    s_col, [], u_col, ...
                                                    C_mat, [], D_mat, b_val);
            start_index = start_index + ny;
            start_index_v = start_index_v + ny;
        end
        if mpc.y_cnstr.max_limit
            %Ineq [C D I -I]*[s u g v]' = y_max-Dd*d
            row = start_index:start_index+ny-1;
            row_v = start_index_v:start_index_v+mpc.ny-1;

            mpc.y_cnstr.max_row_k = row;
            mpc.y_cnstr.max_row_v_k = row_v;

            gi_index_k = mpc.g_index_k(row,k);
            mpc.y_cnstr.g_max_index_k = [mpc.y_cnstr.g_max_index_k gi_index_k];
            mpc.y_cnstr.ineq_max_row = row;

            vi_index_k = mpc.v_index_k(row_v,k);
            mpc.y_cnstr.v_max_index_k = [mpc.y_cnstr.v_max_index_k vi_index_k];
            vi_k = [vi_k; row'];

            C_mat = mpc.C;
            D_mat = mpc.D;
            b_val = mpc.y_cnstr.max;

            [Ai, bi] = appendGeneralizedConstraint(Ai, bi, row,...
                                                    s_col, [], u_col, ...
                                                    C_mat, [], D_mat, b_val);
            start_index = start_index + ny;
            start_index_v = start_index_v + ny;
        end
    end

    if mpc.has_h_cnstr
        if mpc.h_cnstr.min_limit
            %Ineq [-C -Dsu -D I -I]*[s su u g v]' = -h_min+Dd*d
            row = start_index:start_index+nh-1;
            row_v = start_index_v:start_index_v+nh-1;

            mpc.h_cnstr.min_row_k = row;
            mpc.h_cnstr.min_row_v_k = row_v;

            gi_index_k = mpc.g_index_k(row,k);
            mpc.h_cnstr.g_min_index_k = [mpc.h_cnstr.g_min_index_k gi_index_k];
            mpc.h_cnstr.ineq_min_row = row;

            vi_index_k = mpc.v_index_k(row_v,k);
            mpc.h_cnstr.v_min_index_k = [mpc.h_cnstr.v_min_index_k vi_index_k];
            vi_k = [vi_k; row'];

            C_mat = -mpc.Ch;
            Dsu_mat = -mpc.Dsuh;
            D_mat = -mpc.Dh;
            b_val = -mpc.h_cnstr.min;

            [Ai, bi] = appendGeneralizedConstraint(Ai, bi, row,...
                                                    s_col, su_col, u_col, ...
                                                    C_mat, Dsu_mat, D_mat, b_val);
            start_index = start_index + nh;
            start_index_v = start_index_v + nh;
        end
        if mpc.h_cnstr.max_limit
            %Ineq [C Dsu D I -I]*[s su u g v]' = y_max-Dd*d
            row = start_index:start_index+nh-1;
            row_v = start_index_v:start_index_v+nh-1;

            mpc.h_cnstr.max_row_k = row;
            mpc.h_cnstr.max_row_v_k = row_v;

            gi_index_k = mpc.g_index_k(row,k);
            mpc.h_cnstr.g_max_index_k = [mpc.h_cnstr.g_max_index_k gi_index_k];
            mpc.h_cnstr.ineq_max_row = row;

            vi_index_k = mpc.v_index_k(row_v,k);
            mpc.h_cnstr.v_max_index_k = [mpc.h_cnstr.v_max_index_k vi_index_k];
            vi_k = [vi_k; row'];

            C_mat = mpc.Ch;
            Dsu_mat = mpc.Dsuh;
            D_mat = mpc.Dh;
            b_val = mpc.h_cnstr.max;

            [Ai, bi] = appendGeneralizedConstraint(Ai, bi, row,...
                                                    s_col, su_col, u_col, ...
                                                    C_mat, Dsu_mat, D_mat, b_val);
            start_index = start_index + nh;
            start_index_v = start_index_v + nh;
        end
    end

    Ai_k(:,:,k) = Ai;
    bi_k(:,k) = bi;
    mpc.vi_k = vi_k;
end
mpc.Ai_k = Ai_k;
mpc.bi_k = bi_k;

start_index = 1;
Ai_ter = zeros(mpc.ng_k(mpc.N+1),nse);
bi_ter = zeros(mpc.ng_k(mpc.N+1),1);
s_col = 1:mpc.nx;
if mpc.has_s_cnstr
    if mpc.s_cnstr.min_limit
        %Ineq [-I I -I]*[sk gk vk]'=-s_min
        row = start_index:start_index+mpc.nx-1;

        mpc.s_cnstr.min_row_ter = row;

        gi_index_k = mpc.g_index_ter(row);
        mpc.s_cnstr.g_min_index_k = [mpc.s_cnstr.g_min_index_k gi_index_k];

        vi_index_k = mpc.v_index_ter(row);
        mpc.s_cnstr.v_min_index_k = [mpc.s_cnstr.v_min_index_k vi_index_k];

        b_val = -mpc.s_cnstr.min;

        [Ai_ter, bi_ter] = appendGeneralizedConstraint(Ai_ter, bi_ter, row,...
                                                        s_col, [], [], ...
                                                        -eye(nx), [], [], b_val);
        start_index = start_index + nx;
    end
    if mpc.s_cnstr.max_limit
        %Ineq [I I -I]*[sk gk vk]'=s_max
        row = start_index:start_index+mpc.nx-1;

        mpc.s_cnstr.max_row_ter = row;

        gi_index_k = mpc.g_index_ter(row);
        mpc.s_cnstr.g_min_index_k = [mpc.s_cnstr.g_min_index_k gi_index_k];

        vi_index_k = mpc.v_index_ter(row);
        mpc.s_cnstr.v_min_index_k = [mpc.s_cnstr.v_min_index_k vi_index_k];

        b_val = mpc.s_cnstr.max;

        [Ai_ter, bi_ter] = appendGeneralizedConstraint(Ai_ter, bi_ter, row,...
                                                        s_col, [], [], ...
                                                        eye(nx), [], [], b_val);
    end
end

mpc.Ai_ter = Ai_ter;
mpc.bi_ter = bi_ter;


%% dynamics

if mpc.has_du
    mpc.A_kkt = [mpc.A zeros(mpc.nx,mpc.nu);
             zeros(mpc.nu,mpc.nx) zeros(mpc.nu)];
    mpc.B_kkt = [mpc.B;eye(mpc.nu)];
else
    mpc.A_kkt = mpc.A;
    mpc.B_kkt = mpc.B;
end

end



