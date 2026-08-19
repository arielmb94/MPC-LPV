function mpc = genEqualities(mpc)

nx = mpc.nx;
nu = mpc.nu;
nse = mpc.nse;
nvar_k = nse+nu;
mpc.nvar_k = nvar_k;
ny_0 = mpc.ny_0;
ny = mpc.ny;
ny_ter = mpc.ny_ter;
nh_0 = mpc.nh_0;
nh = mpc.nh;
nh_ter = mpc.nh_ter;
N = mpc.N;

s_col = 1:nx;
su_col = nx+1:nse;
u_col = nse+1:nse+nu;

mpc.s_col = s_col;
mpc.su_col = su_col;
mpc.se_col = [s_col su_col];
mpc.u_col = u_col;

mpc.rp_0 = zeros(nse,1);
mpc.rp_k = zeros(nse,mpc.N-1);
mpc.beq_0 = zeros(nx,1);
mpc.beq_k = zeros(nx,mpc.N-1);

mpc.ri_0 = zeros(mpc.ng_k(1),1);
mpc.ri_k = zeros(mpc.ng_k(2),mpc.N-1);
mpc.ri_ter = zeros(mpc.ng_k(3),1);

mpc.ri_hat_0 = zeros(mpc.ng_k(1),1);
mpc.ri_hat_k = zeros(mpc.ng_k(2),mpc.N-1);
mpc.ri_hat_ter = zeros(mpc.ng_k(3),1);

mpc.S_0 = zeros(mpc.ng_k(1),1);
mpc.S_k = zeros(mpc.ng_k(2),mpc.N-1);
mpc.S_ter = zeros(mpc.ng_k(3),1);

mpc.iS_0 = zeros(mpc.ng_k(1),1);
mpc.iS_k = zeros(mpc.ng_k(2),mpc.N-1);
mpc.iS_ter = zeros(mpc.ng_k(3),1);


%% k = 0
v_rows_0 = [];

Ai_0 = zeros(mpc.ng_k(1),nu);
start_index = 1;

u_col = 1:mpc.nu;

if mpc.has_u_cnstr
    if mpc.u_cnstr.min_limit
        % Lower-bound Jacobian block: -I on u.
        row = start_index:start_index+nu-1;

        gi_index_k = mpc.g_index_0(row);
        mpc.u_cnstr.g_min_index_k = gi_index_k;
        mpc.u_cnstr.min_ineqRow_0 = row;

        Ai_0 = appendGeneralizedConstraint(Ai_0, row,...  
                                                 [], [], u_col, ...
                                                 [], [], -eye(nu));

        start_index = start_index + nu;
    end
    if mpc.u_cnstr.max_limit
        % Upper-bound Jacobian block: I on u.
        row = start_index:start_index+nu-1;

        gi_index_k = mpc.g_index_0(row);
        mpc.u_cnstr.g_max_index_k = gi_index_k;
        mpc.u_cnstr.max_ineqRow_0 = row;

        Ai_0 = appendGeneralizedConstraint(Ai_0, row,...  
                                                 [], [], u_col, ...
                                                 [], [], eye(nu));

        start_index = start_index + nu;
    end
end

if mpc.has_du_cnstr
    if mpc.du_cnstr.min_limit
        % Lower-bound Jacobian block: -I on u.
        row = start_index:start_index+nu-1;

        gi_index_k = mpc.g_index_0(row);
        mpc.du_cnstr.g_min_index_k = gi_index_k;
        mpc.du_cnstr.min_ineqRow_0 = row;

        Ai_0 = appendGeneralizedConstraint(Ai_0, row,...  
                                                 [], [], u_col, ...
                                                 [], [], -eye(nu));

        start_index = start_index + nu;
    end
    if mpc.du_cnstr.max_limit
        % Upper-bound Jacobian block: I on u.
        row = start_index:start_index+nu-1;

        gi_index_k = mpc.g_index_0(row);
        mpc.du_cnstr.g_max_index_k = gi_index_k;
        mpc.du_cnstr.max_ineqRow_0 = row;

        Ai_0 = appendGeneralizedConstraint(Ai_0, row,...  
                                                 [], [], u_col, ...
                                                 [], [], eye(nu));

        start_index = start_index + nu;
    end
end

start_index_v = 1;
if mpc.has_y_cnstr
    if mpc.y_cnstr.min_limit && mpc.y_cnstr.use_k0
        % Lower-bound Jacobian block: -D on u.
        row = start_index:start_index+ny_0-1;
        row_v = start_index_v:start_index_v+ny_0-1;

        mpc.y_cnstr.min_ineqRow_0 = row;
        mpc.y_cnstr.min_row_v_0 = row_v;

        gi_index_k = mpc.g_index_0(row);
        mpc.y_cnstr.g_min_index_0 = gi_index_k;

        vi_index_k = mpc.v_index_0(row_v);
        mpc.y_cnstr.v_min_index_0 = vi_index_k;
        v_rows_0 = [v_rows_0; row'];

        D_mat = -mpc.D_0;
        Ai_0 = appendGeneralizedConstraint(Ai_0, row,...  
                                                 [], [], u_col, ...
                                                 [], [], D_mat);
        start_index = start_index + ny_0;
        start_index_v = start_index_v + ny_0;
    end
    if mpc.y_cnstr.max_limit && mpc.y_cnstr.use_k0
        % Upper-bound Jacobian block: D on u.
        row = start_index:start_index+ny_0-1;
        row_v = start_index_v:start_index_v+ny_0-1;

        mpc.y_cnstr.max_ineqRow_0 = row;
        mpc.y_cnstr.max_row_v_0 = row_v;

        gi_index_k = mpc.g_index_0(row);
        mpc.y_cnstr.g_max_index_0 = gi_index_k;

        vi_index_k = mpc.v_index_0(row_v);
        mpc.y_cnstr.v_max_index_0 = vi_index_k;
        v_rows_0 = [v_rows_0; row'];

        D_mat = mpc.D_0;
        Ai_0 = appendGeneralizedConstraint(Ai_0, row,...  
                                                 [], [], u_col, ...
                                                 [], [], D_mat);
        start_index = start_index + ny_0;
        start_index_v = start_index_v + ny_0;
    end
end

if mpc.has_h_cnstr
    if mpc.h_cnstr.min_limit && mpc.h_cnstr.use_k0
        % Lower-bound Jacobian block: -Dh on u.
        row = start_index:start_index+nh_0-1;
        row_v = start_index_v:start_index_v+nh_0-1;

        mpc.h_cnstr.min_ineqRow_0 = row;
        mpc.h_cnstr.min_row_v_0 = row_v;

        gi_index_k = mpc.g_index_0(row);
        mpc.h_cnstr.g_min_index_0 = gi_index_k;

        vi_index_k = mpc.v_index_0(row_v);
        mpc.h_cnstr.v_min_index_0 = vi_index_k;
        v_rows_0 = [v_rows_0; row'];

        D_mat = -mpc.Dh_0;
        Ai_0 = appendGeneralizedConstraint(Ai_0, row,...  
                                                 [], [], u_col, ...
                                                 [], [], D_mat);
        start_index = start_index + nh_0;
        start_index_v = start_index_v + nh_0;
    end
    if mpc.h_cnstr.max_limit && mpc.h_cnstr.use_k0
        % Upper-bound Jacobian block: Dh on u.
        row = start_index:start_index+nh_0-1;
        row_v = start_index_v:start_index_v+nh_0-1;

        mpc.h_cnstr.max_ineqRow_0 = row;
        mpc.h_cnstr.max_row_v_0 = row_v;

        gi_index_k = mpc.g_index_0(row);
        mpc.h_cnstr.g_max_index_0 = gi_index_k;

        vi_index_k = mpc.v_index_0(row_v);
        mpc.h_cnstr.v_max_index_0 = vi_index_k;
        v_rows_0 = [v_rows_0; row'];

        D_mat = mpc.Dh_0;
        Ai_0 = appendGeneralizedConstraint(Ai_0, row,...  
                                                 [], [], u_col, ...
                                                 [], [], D_mat);
        start_index = start_index + nh_0;
        start_index_v = start_index_v + nh_0;
    end
end
mpc.v_rows_0 = v_rows_0;
mpc.Ai_0 = Ai_0;  
%%
s_col = 1:nx;
su_col = nx+1:nse;
u_col = nse+1:nse+nu;

Ai_k = [];
for k = 1:N-1

    Ai = zeros(mpc.ng_k(2),nvar_k);

    v_rows_k = [];
    start_index = 1;
    start_index_v = 1;

    if mpc.has_s_cnstr
        if mpc.s_cnstr.min_limit
            % Lower-bound Jacobian block: -I on s.
            row = start_index:start_index+nx-1;
            row_v = start_index_v:start_index_v+nx-1;

            mpc.s_cnstr.min_ineqRow_k = row;
            mpc.s_cnstr.min_row_v_k = row_v;
           
            gi_index_k = mpc.g_index_k(row,k);
            mpc.s_cnstr.g_min_index_k = [mpc.s_cnstr.g_min_index_k gi_index_k];

            vi_index_k = mpc.v_index_k(row_v,k);
            mpc.s_cnstr.v_min_index_k = [mpc.s_cnstr.v_min_index_k vi_index_k];
            v_rows_k = [v_rows_k; row'];

            Ai = appendGeneralizedConstraint(Ai, row,...  
                                                   s_col, [], [], ...
                                                   -eye(nx), [], []);
            start_index = start_index + nx;
            start_index_v = start_index_v + nx;
        end
        if mpc.s_cnstr.max_limit
            % Upper-bound Jacobian block: I on s.
            row = start_index:start_index+nx-1;
            row_v = start_index_v:start_index_v+nx-1;

            mpc.s_cnstr.max_ineqRow_k = row;
            mpc.s_cnstr.max_row_v_k = row_v;

            gi_index_k = mpc.g_index_k(row,k);
            mpc.s_cnstr.g_max_index_k = [mpc.s_cnstr.g_max_index_k gi_index_k];

            vi_index_k = mpc.v_index_k(row_v,k);
            mpc.s_cnstr.v_max_index_k = [mpc.s_cnstr.v_max_index_k vi_index_k];
            v_rows_k = [v_rows_k; row'];

            Ai = appendGeneralizedConstraint(Ai, row,...  
                                                   s_col, [], [], ...
                                                   eye(nx), [], []);
            start_index = start_index + nx;
            start_index_v = start_index_v + nx;
        end
    end

    if mpc.has_u_cnstr
        if mpc.u_cnstr.min_limit
            % Lower-bound Jacobian block: -I on u.
            row = start_index:start_index+nu-1;

            mpc.u_cnstr.min_ineqRow_k = row;

            gi_index_k = mpc.g_index_k(row,k);
            mpc.u_cnstr.g_min_index_k = [mpc.u_cnstr.g_min_index_k gi_index_k];

            Ai = appendGeneralizedConstraint(Ai, row,...
                [], [], u_col, ...
                [], [], -eye(nu));

            start_index = start_index + nu;
        end
        if mpc.u_cnstr.max_limit
            % Upper-bound Jacobian block: I on u.
            row = start_index:start_index+nu-1;

            mpc.u_cnstr.max_ineqRow_k = row;

            gi_index_k = mpc.g_index_k(row,k);
            mpc.u_cnstr.g_max_index_k = [mpc.u_cnstr.g_max_index_k gi_index_k];

            Ai = appendGeneralizedConstraint(Ai, row,...
                [], [], u_col, ...
                [], [], eye(nu));

            start_index = start_index + nu;
        end
    end

    if mpc.has_du_cnstr
        if mpc.du_cnstr.min_limit
            % Lower-bound Jacobian block: I on su and -I on u.
            row = start_index:start_index+nu-1;

            mpc.du_cnstr.min_ineqRow_k = row;

            gi_index_k = mpc.g_index_k(row,k);
            mpc.du_cnstr.g_min_index_k = [mpc.du_cnstr.g_min_index_k gi_index_k];

            Ai = appendGeneralizedConstraint(Ai, row,...
                [], su_col, u_col, ...
                [], eye(nu), -eye(nu));

            start_index = start_index + nu;
        end
        if mpc.du_cnstr.max_limit
            % Upper-bound Jacobian block: -I on su and I on u.
            row = start_index:start_index+nu-1;

            mpc.du_cnstr.max_ineqRow_k = row;

            gi_index_k = mpc.g_index_k(row,k);
            mpc.du_cnstr.g_max_index_k = [mpc.du_cnstr.g_max_index_k gi_index_k];

            Ai = appendGeneralizedConstraint(Ai, row,...
                [], su_col, u_col, ...
                [], -eye(nu), eye(nu));

            start_index = start_index + nu;
        end
    end

    if mpc.has_y_cnstr
        if mpc.y_cnstr.min_limit
            % Lower-bound Jacobian block: -C on s and -D on u.
            row = start_index:start_index+ny-1;
            row_v = start_index_v:start_index_v+ny-1;

            mpc.y_cnstr.min_ineqRow_k = row;
            mpc.y_cnstr.min_row_v_k = row_v;

            gi_index_k = mpc.g_index_k(row,k);
            mpc.y_cnstr.g_min_index_k = [mpc.y_cnstr.g_min_index_k gi_index_k];

            vi_index_k = mpc.v_index_k(row_v,k);
            mpc.y_cnstr.v_min_index_k = [mpc.y_cnstr.v_min_index_k vi_index_k];
            v_rows_k = [v_rows_k; row'];

            C_mat = [];
            D_mat = [];

            if mpc.y_cnstr.use_s, C_mat = -mpc.C(:,:,k); end
            if mpc.y_cnstr.use_u, D_mat = -mpc.D(:,:,k); end
            Ai = appendGeneralizedConstraint(Ai, row,...
                                                    s_col, [], u_col, ...
                                                    C_mat, [], D_mat);
            start_index = start_index + ny;
            start_index_v = start_index_v + ny;
        end
        if mpc.y_cnstr.max_limit
            % Upper-bound Jacobian block: C on s and D on u.
            row = start_index:start_index+ny-1;
            row_v = start_index_v:start_index_v+ny-1;

            mpc.y_cnstr.max_ineqRow_k = row;
            mpc.y_cnstr.max_row_v_k = row_v;

            gi_index_k = mpc.g_index_k(row,k);
            mpc.y_cnstr.g_max_index_k = [mpc.y_cnstr.g_max_index_k gi_index_k];

            vi_index_k = mpc.v_index_k(row_v,k);
            mpc.y_cnstr.v_max_index_k = [mpc.y_cnstr.v_max_index_k vi_index_k];
            v_rows_k = [v_rows_k; row'];

            C_mat = [];
            D_mat = [];

            if mpc.y_cnstr.use_s, C_mat = mpc.C(:,:,k); end
            if mpc.y_cnstr.use_u, D_mat = mpc.D(:,:,k); end
            Ai = appendGeneralizedConstraint(Ai, row,...
                                                    s_col, [], u_col, ...
                                                    C_mat, [], D_mat);
            start_index = start_index + ny;
            start_index_v = start_index_v + ny;
        end
    end

    if mpc.has_h_cnstr
        if mpc.h_cnstr.min_limit
            % Lower-bound Jacobian block: -C on s, -Dsu on su, and -D on u.
            row = start_index:start_index+nh-1;
            row_v = start_index_v:start_index_v+nh-1;

            mpc.h_cnstr.min_ineqRow_k = row;
            mpc.h_cnstr.min_row_v_k = row_v;

            gi_index_k = mpc.g_index_k(row,k);
            mpc.h_cnstr.g_min_index_k = [mpc.h_cnstr.g_min_index_k gi_index_k];

            vi_index_k = mpc.v_index_k(row_v,k);
            mpc.h_cnstr.v_min_index_k = [mpc.h_cnstr.v_min_index_k vi_index_k];
            v_rows_k = [v_rows_k; row'];

            C_mat = [];
            Dsu_mat = [];
            D_mat = [];
            
            if mpc.h_cnstr.use_s, C_mat = -mpc.Ch(:,:,k); end
            if mpc.h_cnstr.use_su, Dsu_mat = -mpc.Dsuh(:,:,k); end
            if mpc.h_cnstr.use_u, D_mat = -mpc.Dh(:,:,k); end
            Ai = appendGeneralizedConstraint(Ai, row,...
                                                    s_col, su_col, u_col, ...
                                                    C_mat, Dsu_mat, D_mat);
            start_index = start_index + nh;
            start_index_v = start_index_v + nh;
        end
        if mpc.h_cnstr.max_limit
            % Upper-bound Jacobian block: C on s, Dsu on su, and D on u.
            row = start_index:start_index+nh-1;
            row_v = start_index_v:start_index_v+nh-1;

            mpc.h_cnstr.max_ineqRow_k = row;
            mpc.h_cnstr.max_row_v_k = row_v;

            gi_index_k = mpc.g_index_k(row,k);
            mpc.h_cnstr.g_max_index_k = [mpc.h_cnstr.g_max_index_k gi_index_k];

            vi_index_k = mpc.v_index_k(row_v,k);
            mpc.h_cnstr.v_max_index_k = [mpc.h_cnstr.v_max_index_k vi_index_k];
            v_rows_k = [v_rows_k; row'];

            C_mat = [];
            Dsu_mat = [];
            D_mat = [];
            
            if mpc.h_cnstr.use_s, C_mat = mpc.Ch(:,:,k); end
            if mpc.h_cnstr.use_su, Dsu_mat = mpc.Dsuh(:,:,k); end
            if mpc.h_cnstr.use_u, D_mat = mpc.Dh(:,:,k); end
            Ai = appendGeneralizedConstraint(Ai, row,...
                                                    s_col, su_col, u_col, ...
                                                    C_mat, Dsu_mat, D_mat);
            start_index = start_index + nh;
            start_index_v = start_index_v + nh;
        end
    end

    Ai_k(:,:,k) = Ai;
    mpc.v_rows_k = v_rows_k;
end
mpc.Ai_k = Ai_k;
%%
start_index = 1;
Ai_ter = zeros(mpc.ng_k(3),nse);

if mpc.has_s_cnstr
    if mpc.s_cnstr.min_limit
        % Lower-bound Jacobian block: -I on s.
        row = start_index:start_index+nx-1;

        mpc.s_cnstr.min_ineqRow_ter = row;

        gi_index_k = mpc.g_index_ter(row);
        mpc.s_cnstr.g_min_index_k = [mpc.s_cnstr.g_min_index_k gi_index_k];

        vi_index_k = mpc.v_index_ter(row);
        mpc.s_cnstr.v_min_index_k = [mpc.s_cnstr.v_min_index_k vi_index_k];

        Ai_ter = appendGeneralizedConstraint(Ai_ter, row,...
                                                        s_col, [], [], ...
                                                        -eye(nx), [], []);
        start_index = start_index + nx;
    end
    if mpc.s_cnstr.max_limit
        % Upper-bound Jacobian block: I on s.
        row = start_index:start_index+nx-1;

        mpc.s_cnstr.max_ineqRow_ter = row;

        gi_index_k = mpc.g_index_ter(row);
        mpc.s_cnstr.g_min_index_k = [mpc.s_cnstr.g_min_index_k gi_index_k];

        vi_index_k = mpc.v_index_ter(row);
        mpc.s_cnstr.v_min_index_k = [mpc.s_cnstr.v_min_index_k vi_index_k];

        Ai_ter = appendGeneralizedConstraint(Ai_ter, row,...
                                                        s_col, [], [], ...
                                                        eye(nx), [], []);
        start_index = start_index + nx;
    end
end

if mpc.has_y_cnstr
    if mpc.y_cnstr.min_limit && mpc.y_cnstr.use_ter
        % Lower-bound Jacobian block: -C on s.
        row = start_index:start_index+ny_ter-1;

        mpc.y_cnstr.min_ineqRow_ter = row;

        gi_index_k = mpc.g_index_ter(row);
        mpc.y_cnstr.g_min_index_ter = gi_index_k;

        vi_index_k = mpc.v_index_ter(row);
        mpc.y_cnstr.v_min_index_ter = vi_index_k;

        C_mat = -mpc.C_ter;
        Ai_ter = appendGeneralizedConstraint(Ai_ter, row,...
                                               s_col, [], [], ...
                                               C_mat, [], []);
        
        start_index = start_index + ny_ter;
    end
    if mpc.y_cnstr.max_limit && mpc.y_cnstr.use_ter
        % Upper-bound Jacobian block: C on s.
        row = start_index:start_index+ny_ter-1;

        mpc.y_cnstr.max_ineqRow_ter = row;

        gi_index_k = mpc.g_index_ter(row);
        mpc.y_cnstr.g_max_index_ter = gi_index_k;

        vi_index_k = mpc.v_index_ter(row);
        mpc.y_cnstr.v_max_index_ter = vi_index_k;

        C_mat = mpc.C_ter;
        Ai_ter = appendGeneralizedConstraint(Ai_ter, row,...
                                               s_col, [], [], ...
                                               C_mat, [], []);
        start_index = start_index + ny_ter;
    end
end

if mpc.has_h_cnstr
    if mpc.h_cnstr.min_limit && mpc.h_cnstr.use_ter
        % Lower-bound Jacobian block: -Ch on s.
        row = start_index:start_index+nh_ter-1;

        mpc.h_cnstr.min_ineqRow_ter = row;

        gi_index_k = mpc.g_index_ter(row);
        mpc.h_cnstr.g_min_index_ter = gi_index_k;

        vi_index_k = mpc.v_index_ter(row);
        mpc.h_cnstr.v_min_index_ter = vi_index_k;

        C_mat = -mpc.Ch_ter;
        Ai_ter = appendGeneralizedConstraint(Ai_ter, row,...
                                               s_col, [], [], ...
                                               C_mat, [], []);
        
        start_index = start_index + nh_ter;
    end
    if mpc.h_cnstr.max_limit && mpc.h_cnstr.use_ter
        % Upper-bound Jacobian block: Ch on s.
        row = start_index:start_index+nh_ter-1;

        mpc.h_cnstr.max_ineqRow_ter = row;

        gi_index_k = mpc.g_index_ter(row);
        mpc.h_cnstr.g_max_index_ter = gi_index_k;

        vi_index_k = mpc.v_index_ter(row);
        mpc.h_cnstr.v_max_index_ter = vi_index_k;

        C_mat = mpc.Ch_ter;
        Ai_ter = appendGeneralizedConstraint(Ai_ter, row,...
                                               s_col, [], [], ...
                                               C_mat, [], []);
        start_index = start_index + nh_ter;
    end
end

mpc.Ai_ter = Ai_ter;


%% dynamics

if mpc.has_du
    mpc.A_kkt = zeros(mpc.nse,mpc.nse,mpc.N-1);
    mpc.A_kkt(mpc.s_col,mpc.s_col,:) = mpc.A(:,:,2:mpc.N);

    mpc.B_kkt_0 = [mpc.B(:,:,1);eye(mpc.nu)];
    mpc.B_kkt = zeros(mpc.nse,mpc.nu,mpc.N-1);
    mpc.B_kkt(mpc.s_col,:,:) = mpc.B(:,:,2:mpc.N);
    mpc.B_kkt(mpc.su_col,:,:) = eye(mpc.nu).*ones(1,1,N-1);

else
    mpc.A_kkt = mpc.A(:,:,2:mpc.N);

    mpc.B_kkt_0 = mpc.B(:,:,1);
    mpc.B_kkt = mpc.B(:,:,2:mpc.N);
end

end
