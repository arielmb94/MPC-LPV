function mpc = build_index(mpc)

if mpc.has_s_cnstr
ns_min = mpc.s_cnstr.min_limit*mpc.nx;
ns_max = mpc.s_cnstr.max_limit*mpc.nx;
else
ns_min = 0;
ns_max = 0;
end

if mpc.has_u_cnstr
nu_min = mpc.u_cnstr.min_limit*mpc.nu;
nu_max = mpc.u_cnstr.max_limit*mpc.nu;
else
nu_min = 0;
nu_max = 0;
end

du = mpc.has_du;
if mpc.has_du_cnstr
ndu_min = mpc.du_cnstr.min_limit*mpc.nu;
ndu_max = mpc.du_cnstr.max_limit*mpc.nu;
else
ndu_min = 0;
ndu_max = 0;
end

if mpc.has_y_cnstr

ny_min_0 = mpc.y_cnstr.min_limit*mpc.ny_0;
ny_max_0 = mpc.y_cnstr.max_limit*mpc.ny_0;    

ny_min = mpc.y_cnstr.min_limit*mpc.ny;
ny_max = mpc.y_cnstr.max_limit*mpc.ny;

ny_min_ter = mpc.y_cnstr.min_limit*mpc.ny_ter;
ny_max_ter = mpc.y_cnstr.max_limit*mpc.ny_ter;

else
ny_min_0 = 0;
ny_max_0 = 0;
ny_min = 0;
ny_max = 0;
ny_min_ter = 0;
ny_max_ter = 0;
end

if mpc.has_h_cnstr
nh_min_0 = mpc.h_cnstr.min_limit*mpc.nh_0;
nh_max_0 = mpc.h_cnstr.max_limit*mpc.nh_0;

nh_min = mpc.h_cnstr.min_limit*mpc.nh;
nh_max = mpc.h_cnstr.max_limit*mpc.nh;

nh_min_ter = mpc.h_cnstr.min_limit*mpc.nh_ter;
nh_max_ter = mpc.h_cnstr.max_limit*mpc.nh_ter;

else
nh_min_0 = 0;
nh_max_0 = 0;
nh_min = 0;
nh_max = 0;
nh_min_ter = 0;
nh_max_ter = 0;
end


%% Optimization variables x
mpc.n = mpc.Nx+mpc.Nu+mpc.Nu*du+...
        (ns_min*mpc.N+ns_max*mpc.N)*2+...
        nu_min*mpc.N+nu_max*mpc.N+...
        (ndu_min*mpc.N+ndu_max*mpc.N)*du+...
        (ny_min*(mpc.N-1)+ny_max*(mpc.N-1))*2+...
         (ny_min_0+ny_max_0+ny_min_ter+ny_max_ter)*2+...
        (nh_min*(mpc.N-1)+nh_max*(mpc.N-1))*2+...
        (nh_min_0+nh_max_0+nh_min_ter+nh_max_ter)*2;

s_index = [];
su_index = [];
u_index = [];
g_index = [];
v_index = [];

s_index_k = [];
su_index_k = [];
u_index_k = [];

g_index_0 = [];
g_index_k = [];
g_index_ter = [];
v_index_0 = [];
v_index_k = [];
v_index_ter = [];

gs_min_index_k = [];
gs_max_index_k = [];
gu_min_index_k = [];
gu_max_index_k = [];
gdu_min_index_k = [];
gdu_max_index_k = [];
gy_min_index_k = [];
gy_max_index_k = [];
gh_min_index_k = [];
gh_max_index_k = [];

vs_min_index_k = [];
vs_max_index_k = [];
vy_min_index_k = [];
vy_max_index_k = [];
vh_min_index_k = [];
vh_max_index_k = [];

%% x
start_index = 1;

% k = 0
dim = mpc.nu;   % u

% u
[start_index,u_index,u_index_k] = expand_index(dim,start_index,1,u_index);

for k = 1:mpc.N-1

dim = [mpc.nx mpc.nu*du mpc.nu];                                                  % s su u

% s
[start_index,s_index,s_index_k_vec] = expand_index(dim,start_index,1,s_index);

% su
[start_index,su_index,su_index_k_vec] = expand_index(dim,start_index,2,su_index);

% u
[start_index,u_index,u_index_k_vec] = expand_index(dim,start_index,3,u_index);

s_index_k = [s_index_k s_index_k_vec];
su_index_k = [su_index_k su_index_k_vec];
u_index_k = [u_index_k u_index_k_vec];

end

% k = N
dim = [mpc.nx mpc.nu*du];           % s su
% s
[start_index,s_index,s_index_k_vec] = expand_index(dim,start_index,1,s_index);

% su
[start_index,su_index,su_index_k_vec] = expand_index(dim,start_index,2,su_index);

s_index_k = [s_index_k s_index_k_vec];
su_index_k = [su_index_k su_index_k_vec];

mpc.u_index = u_index;
mpc.s_index = s_index;
mpc.su_index = su_index;

mpc.u_index_k = u_index_k;
mpc.s_index_k = s_index_k;
mpc.su_index_k = su_index_k;
mpc.nse = mpc.nx+du*mpc.nu;
mpc.se_index_k = [s_index_k;su_index_k];

mpc.delta_u = mpc.u_index_k*0;
mpc.delta_se = mpc.se_index_k*0;

mpc.nvar = start_index-1;
mpc.variables_index = [1:mpc.nvar]';

mpc.ru_0 = zeros(mpc.nu,1);
mpc.ru_k = zeros(mpc.nu,mpc.N-1);
mpc.rse_k = zeros(mpc.nse,mpc.N-1);
mpc.rse_ter = zeros(mpc.nse,1);

mpc.ru_hat_0 = zeros(mpc.nu,1);
mpc.ru_hat_k = zeros(mpc.nu,mpc.N-1);
mpc.rse_hat_k = zeros(mpc.nse,mpc.N-1);
mpc.rse_hat_ter = zeros(mpc.nse,1);

%% g

% k = 0
dim = [nu_min nu_max ndu_min ndu_max ny_min_0 ny_max_0 nh_min_0 nh_max_0];  % gu gdu gy_0 gh_0

% gu
[start_index,g_index,gu_min_index_k] = expand_index(dim,start_index,1,g_index);
[start_index,g_index,gu_max_index_k] = expand_index(dim,start_index,2,g_index);

% gdu
[start_index,g_index,gdu_min_index_k] = expand_index(dim,start_index,3,g_index);
[start_index,g_index,gdu_max_index_k] = expand_index(dim,start_index,4,g_index);

% gy_0
[start_index,g_index,gy_min_index_k] = expand_index(dim,start_index,5,g_index);
[start_index,g_index,gy_max_index_k] = expand_index(dim,start_index,6,g_index);

% gh_0
[start_index,g_index,gh_min_index_k] = expand_index(dim,start_index,7,g_index);
[start_index,g_index,gh_max_index_k] = expand_index(dim,start_index,8,g_index);

g_index_0 = [gu_min_index_k;gu_max_index_k;
            gdu_min_index_k;gdu_max_index_k;
            gy_min_index_k;gy_max_index_k;
            gh_min_index_k;gh_max_index_k];

for k = 1:mpc.N-1

dim = [ns_min ns_max nu_min nu_max ndu_min ndu_max ny_min ny_max nh_min nh_max];   % gs gu gdu gy gh

% gs
[start_index,g_index,gs_min_index_k] = expand_index(dim,start_index,1,g_index);
[start_index,g_index,gs_max_index_k] = expand_index(dim,start_index,2,g_index);

% gu
[start_index,g_index,gu_min_index_k] = expand_index(dim,start_index,3,g_index);
[start_index,g_index,gu_max_index_k] = expand_index(dim,start_index,4,g_index);

% gdu
[start_index,g_index,gdu_min_index_k] = expand_index(dim,start_index,5,g_index);
[start_index,g_index,gdu_max_index_k] = expand_index(dim,start_index,6,g_index);

% gy
[start_index,g_index,gy_min_index_k] = expand_index(dim,start_index,7,g_index);
[start_index,g_index,gy_max_index_k] = expand_index(dim,start_index,8,g_index);

% gh
[start_index,g_index,gh_min_index_k] = expand_index(dim,start_index,9,g_index);
[start_index,g_index,gh_max_index_k] = expand_index(dim,start_index,10,g_index);

g_index_k_vec = [gs_min_index_k;gs_max_index_k;
                 gu_min_index_k;gu_max_index_k;
                 gdu_min_index_k;gdu_max_index_k;
                 gy_min_index_k;gy_max_index_k;
                 gh_min_index_k;gh_max_index_k];

g_index_k = [g_index_k g_index_k_vec];

end

% k = N
dim = [ns_min ns_max ny_min_ter ny_max_ter nh_min_ter nh_max_ter];  % gs gy_ter gh_ter

% gs
[start_index,g_index,gs_min_index_k] = expand_index(dim,start_index,1,g_index);
[start_index,g_index,gs_max_index_k] = expand_index(dim,start_index,2,g_index);

% gy
[start_index,g_index,gy_min_index_k] = expand_index(dim,start_index,3,g_index);
[start_index,g_index,gy_max_index_k] = expand_index(dim,start_index,4,g_index);

% gh
[start_index,g_index,gh_min_index_k] = expand_index(dim,start_index,5,g_index);
[start_index,g_index,gh_max_index_k] = expand_index(dim,start_index,6,g_index);

g_index_ter = [gs_min_index_k;gs_max_index_k;
               gy_min_index_k;gy_max_index_k;
               gh_min_index_k;gh_max_index_k];

mpc.g_index = g_index;

mpc.g_index_0 = g_index_0;
mpc.g_index_k = g_index_k;
mpc.g_index_ter = g_index_ter;

mpc.g_0 = g_index_0*0;
mpc.g_k = g_index_k*0;
mpc.g_ter = g_index_ter*0;

mpc.rg_0 = g_index_0*0;
mpc.rg_k = g_index_k*0;
mpc.rg_ter = g_index_ter*0;

mpc.mu_i_0 = g_index_0*0;
mpc.mu_i_k = g_index_k*0;
mpc.mu_i_ter = g_index_ter*0;

mpc.g2_0 = g_index_0*0;
mpc.g2_k = g_index_k*0;
mpc.g2_ter = g_index_ter*0;

mpc.delta_g_0 = g_index_0*0;
mpc.delta_g_k = g_index_k*0;
mpc.delta_g_ter = g_index_ter*0;

%% v

% k = 0
dim = [ny_min_0 ny_max_0 nh_min_0 nh_max_0];  % vy_0 vh_0

% vy_0
[start_index,v_index,vy_min_index_k] = expand_index(dim,start_index,1,v_index);
[start_index,v_index,vy_max_index_k] = expand_index(dim,start_index,2,v_index);

% vh_0
[start_index,v_index,vh_min_index_k] = expand_index(dim,start_index,3,v_index);
[start_index,v_index,vh_max_index_k] = expand_index(dim,start_index,4,v_index);

v_index_0 = [vy_min_index_k;vy_max_index_k;
             vh_min_index_k;vh_max_index_k];

for k = 1:mpc.N-1

dim = [ns_min ns_max ny_min ny_max nh_min nh_max];  % vs vy vh

% vs
[start_index,v_index,vs_min_index_k] = expand_index(dim,start_index,1,v_index);
[start_index,v_index,vs_max_index_k] = expand_index(dim,start_index,2,v_index);

% vy
[start_index,v_index,vy_min_index_k] = expand_index(dim,start_index,3,v_index);
[start_index,v_index,vy_max_index_k] = expand_index(dim,start_index,4,v_index);

% vh
[start_index,v_index,vh_min_index_k] = expand_index(dim,start_index,5,v_index);
[start_index,v_index,vh_max_index_k] = expand_index(dim,start_index,6,v_index);

v_index_k_vec = [vs_min_index_k;vs_max_index_k;
                 vy_min_index_k;vy_max_index_k;
                 vh_min_index_k;vh_max_index_k];

v_index_k = [v_index_k v_index_k_vec];

end

% k = N
dim = [ns_min ns_max ny_min_ter ny_max_ter nh_min_ter nh_max_ter];         % vs vh_ter

% vs
[start_index,v_index,vs_min_index_k] = expand_index(dim,start_index,1,v_index);
[start_index,v_index,vs_max_index_k] = expand_index(dim,start_index,2,v_index);

% vy
[start_index,v_index,vy_min_index_k] = expand_index(dim,start_index,3,v_index);
[start_index,v_index,vy_max_index_k] = expand_index(dim,start_index,4,v_index);

% vh
[start_index,v_index,vh_min_index_k] = expand_index(dim,start_index,5,v_index);
[start_index,v_index,vh_max_index_k] = expand_index(dim,start_index,6,v_index);

v_index_ter = [vs_min_index_k;vs_max_index_k;
               vy_min_index_k;vy_max_index_k
               vh_min_index_k;vh_max_index_k];

mpc.v_index = v_index;

mpc.v_index_0 = v_index_0;
mpc.v_index_k = v_index_k;
mpc.v_index_ter = v_index_ter;

mpc.v_0 = v_index_0*0;
mpc.v_k = v_index_k*0;
mpc.v_ter = v_index_ter*0;

mpc.rv_v2_0 = v_index_0*0;
mpc.rv_v2_k = v_index_k*0;
mpc.rv_v2_ter = v_index_ter*0;

mpc.v2_0 = v_index_0*0;
mpc.v2_k = v_index_k*0;
mpc.v2_ter = v_index_ter*0;

mpc.delta_v_0 = v_index_0*0;
mpc.delta_v_k = v_index_k*0;
mpc.delta_v_ter = v_index_ter*0;

end
