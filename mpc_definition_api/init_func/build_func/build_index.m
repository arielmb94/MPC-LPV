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
ny_min = mpc.y_cnstr.min_limit*mpc.ny;
ny_max = mpc.y_cnstr.max_limit*mpc.ny;
else
ny_min = 0;
ny_max = 0;
end

if mpc.has_h_cnstr
nh_min = mpc.h_cnstr.min_limit*mpc.nh;
nh_max = mpc.h_cnstr.max_limit*mpc.nh;

nh_min_0 = nh_min*any(mpc.Dh);
nh_max_0 = nh_max*any(mpc.Dh);
else
nh_min = 0;
nh_max = 0;
nh_min_0 = 0;
nh_max_0 = 0;
end


%% Optimization variables x
mpc.n = mpc.Nx+mpc.Nu+(mpc.Nu-mpc.nu)*du+...
        (ns_min*mpc.N+ns_max*mpc.N)*2+...
        nu_min*mpc.N+nu_max*mpc.N+...
        (ndu_min*mpc.N+ndu_max*mpc.N)*du+...
        (ny_min*(mpc.N-1)+ny_max*(mpc.N-1))*2+...
        (nh_min*(mpc.N-1)+nh_max*(mpc.N-1))*2+...
        nh_min_0*2+nh_max_0*2;

x = zeros(mpc.n,1);

s_index = [];
su_index = [];
u_index = [];
g_index = [];
v_index = [];

s_index_k = [];
su_index_k = [];
u_index_k = [];

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

start_index = 1;

% k = 0
dim = [mpc.nu ...                                           % u
       nu_min nu_max ndu_min ndu_max nh_min_0 nh_max_0,...  % gu gdu gh_0
       nh_min_0 nh_max_0];                                  % vh_0
% u
[start_index,u_index,u_index_k] = expand_index(dim,start_index,1,u_index,u_index_k);
% gu
[start_index,g_index,gu_min_index_k] = expand_index(dim,start_index,2,g_index,gu_min_index_k);
[start_index,g_index,gu_max_index_k] = expand_index(dim,start_index,3,g_index,gu_max_index_k);
% gdu
[start_index,g_index,gdu_min_index_k] = expand_index(dim,start_index,4,g_index,gdu_min_index_k);
[start_index,g_index,gdu_max_index_k] = expand_index(dim,start_index,5,g_index,gdu_max_index_k);
% gh_0
[start_index,g_index,gh_min_index_k] = expand_index(dim,start_index,6,g_index,gh_min_index_k);
[start_index,g_index,gh_max_index_k] = expand_index(dim,start_index,7,g_index,gh_max_index_k);
% vh_0
[start_index,v_index,vh_min_index_k] = expand_index(dim,start_index,8,v_index,vh_min_index_k);
[start_index,v_index,vh_max_index_k] = expand_index(dim,start_index,9,v_index,vh_max_index_k);

for k = 1:mpc.N-1

dim = [mpc.nx mpc.nu*du mpc.nu ...                                                  % s su u
       ns_min ns_max nu_min nu_max ndu_min ndu_max ny_min ny_max nh_min nh_max...   % gs gu gdu gy gh
       ns_min ns_max ny_min ny_max nh_min nh_max];                                  % vs vy vh

% s
[start_index,s_index,s_index_k] = expand_index(dim,start_index,1,s_index,s_index_k);

% su
[start_index,su_index,su_index_k] = expand_index(dim,start_index,2,su_index,su_index_k);

% u
[start_index,u_index,u_index_k] = expand_index(dim,start_index,3,u_index,u_index_k);

% gs
[start_index,g_index,gs_min_index_k] = expand_index(dim,start_index,4,g_index,gs_min_index_k);
[start_index,g_index,gs_max_index_k] = expand_index(dim,start_index,5,g_index,gs_max_index_k);

% gu
[start_index,g_index,gu_min_index_k] = expand_index(dim,start_index,6,g_index,gu_min_index_k);
[start_index,g_index,gu_max_index_k] = expand_index(dim,start_index,7,g_index,gu_max_index_k);

% gdu
[start_index,g_index,gdu_min_index_k] = expand_index(dim,start_index,8,g_index,gdu_min_index_k);
[start_index,g_index,gdu_max_index_k] = expand_index(dim,start_index,9,g_index,gdu_max_index_k);

% gy
[start_index,g_index,gy_min_index_k] = expand_index(dim,start_index,10,g_index,gy_min_index_k);
[start_index,g_index,gy_max_index_k] = expand_index(dim,start_index,11,g_index,gy_max_index_k);

% gh
[start_index,g_index,gh_min_index_k] = expand_index(dim,start_index,12,g_index,gh_min_index_k);
[start_index,g_index,gh_max_index_k] = expand_index(dim,start_index,13,g_index,gh_max_index_k);

% vs
[start_index,v_index,vs_min_index_k] = expand_index(dim,start_index,14,v_index,vs_min_index_k);
[start_index,v_index,vs_max_index_k] = expand_index(dim,start_index,15,v_index,vs_max_index_k);

% vy
[start_index,v_index,vy_min_index_k] = expand_index(dim,start_index,16,v_index,vy_min_index_k);
[start_index,v_index,vy_max_index_k] = expand_index(dim,start_index,17,v_index,vy_max_index_k);

% vh
[start_index,v_index,vh_min_index_k] = expand_index(dim,start_index,18,v_index,vh_min_index_k);
[start_index,v_index,vh_max_index_k] = expand_index(dim,start_index,19,v_index,vh_max_index_k);

end

% k = N
dim = [mpc.nx ...           % s
    ns_min ns_max ...       % gs
    ns_min ns_max];         % vs

% s
[start_index,s_index,s_index_k] = expand_index(dim,start_index,1,s_index,s_index_k);

% gs
[start_index,g_index,gs_min_index_k] = expand_index(dim,start_index,2,g_index,gs_min_index_k);
[start_index,g_index,gs_max_index_k] = expand_index(dim,start_index,3,g_index,gs_max_index_k);

% vs
[start_index,v_index,vs_min_index_k] = expand_index(dim,start_index,4,v_index,vs_min_index_k);
[start_index,v_index,vs_max_index_k] = expand_index(dim,start_index,5,v_index,vs_max_index_k);


mpc.s_index = s_index';
mpc.s_index_k = [zeros(mpc.nx,1) s_index_k'];
mpc.su_index = su_index';
mpc.su_index_k = [zeros(mpc.nu,1) su_index_k' zeros(mpc.nu,1)];
mpc.u_index = u_index';
mpc.u_index_k = [u_index_k' zeros(mpc.nu,1)];
mpc.g_index = g_index';
mpc.v_index = v_index';

x(mpc.g_index) = 1;
x(mpc.v_index) = 1;
mpc.slack_index = find(x==1);
mpc.variables_index = find(x==0);

mpc.g = x(mpc.g_index);
mpc.v = x(mpc.v_index);
mpc.slacks = x(mpc.slack_index);

if mpc.has_s_cnstr
    if mpc.s_cnstr.min_limit
    mpc.s_cnstr.g_min_index_k = [zeros(mpc.nx,1) gs_min_index_k'];
    mpc.s_cnstr.v_min_index_k = [zeros(mpc.nx,1) vs_min_index_k'];
    end
    if mpc.s_cnstr.max_limit
    mpc.s_cnstr.g_max_index_k = [zeros(mpc.nx,1) gs_max_index_k'];
    mpc.s_cnstr.v_max_index_k = [zeros(mpc.nx,1) vs_max_index_k'];
    end
end

if mpc.has_u_cnstr
    if mpc.u_cnstr.min_limit
    mpc.u_cnstr.g_min_index_k = [gu_min_index_k' zeros(mpc.nu,1)];
    end
    if mpc.u_cnstr.max_limit
    mpc.u_cnstr.g_max_index_k = [gu_max_index_k' zeros(mpc.nu,1)];
    end
end

if mpc.has_du_cnstr
    if mpc.du_cnstr.min_limit
    mpc.du_cnstr.g_min_index_k = [gdu_min_index_k' zeros(mpc.nu,1)];
    end
    if mpc.du_cnstr.max_limit
    mpc.du_cnstr.g_max_index_k = [gdu_max_index_k' zeros(mpc.nu,1)];
    end
end

if mpc.has_y_cnstr
    if mpc.y_cnstr.min_limit
    mpc.y_cnstr.g_min_index_k = [zeros(mpc.ny,1) gy_min_index_k' zeros(mpc.ny,1)];
    mpc.y_cnstr.v_min_index_k = [zeros(mpc.ny,1) vy_min_index_k' zeros(mpc.ny,1)];
    end
    if mpc.y_cnstr.max_limit
    mpc.y_cnstr.g_max_index_k = [zeros(mpc.ny,1) gy_max_index_k' zeros(mpc.ny,1)];
    mpc.y_cnstr.v_max_index_k = [zeros(mpc.ny,1) vy_max_index_k' zeros(mpc.ny,1)];
    end
end

if mpc.has_h_cnstr
    if mpc.h_cnstr.min_limit
    mpc.h_cnstr.g_min_index_k = [zeros(mpc.nh-nh_min_0,1) gh_min_index_k' zeros(mpc.nh,1)];
    mpc.h_cnstr.v_min_index_k = [zeros(mpc.nh-nh_min_0,1) vh_min_index_k' zeros(mpc.nh,1)];
    end
    if mpc.h_cnstr.max_limit
    mpc.h_cnstr.g_max_index_k = [zeros(mpc.nh-nh_max_0,1) gh_max_index_k' zeros(mpc.nh,1)];
    mpc.h_cnstr.v_max_index_k = [zeros(mpc.nh-nh_max_0,1) vh_max_index_k' zeros(mpc.nh,1)];
    end
end

%% Equality rows 
start_index = 1;

dyn_k = [];
su_dyn_k = [];
s_min_eq_index_k = [];
s_max_eq_index_k = [];
u_min_eq_index_k = [];
u_max_eq_index_k = [];
du_min_eq_index_k = [];
du_max_eq_index_k = [];
y_min_eq_index_k = [];
y_max_eq_index_k = [];
h_min_eq_index_k = [];
h_max_eq_index_k = [];


% k = 0
% Dynamics
[start_index,dyn_k] = expand_equal_index(mpc.nx,start_index,dyn_k);
[start_index,su_dyn_k] = expand_equal_index(mpc.nu*du,start_index,su_dyn_k);

% U constr
[start_index,u_min_eq_index_k] = expand_equal_index(nu_min,start_index,u_min_eq_index_k);
[start_index,u_max_eq_index_k] = expand_equal_index(nu_max,start_index,u_max_eq_index_k);

% DU constr
[start_index,du_min_eq_index_k] = expand_equal_index(ndu_min,start_index,du_min_eq_index_k);
[start_index,du_max_eq_index_k] = expand_equal_index(ndu_max,start_index,du_max_eq_index_k);

% H constr
[start_index,h_min_eq_index_k] = expand_equal_index(nh_min_0,start_index,h_min_eq_index_k);
[start_index,h_max_eq_index_k] = expand_equal_index(nh_min_0,start_index,h_max_eq_index_k);

for k = 1:mpc.N-1

    % Dynamics
    [start_index,dyn_k] = expand_equal_index(mpc.nx,start_index,dyn_k);
    if k ~= mpc.N-1
        [start_index,su_dyn_k] = expand_equal_index(mpc.nu*du,start_index,su_dyn_k);
    end

    % S constr
    [start_index,s_min_eq_index_k] = expand_equal_index(ns_min,start_index,s_min_eq_index_k);
    [start_index,s_max_eq_index_k] = expand_equal_index(ns_max,start_index,s_max_eq_index_k);

    % U constr
    [start_index,u_min_eq_index_k] = expand_equal_index(nu_min,start_index,u_min_eq_index_k);
    [start_index,u_max_eq_index_k] = expand_equal_index(nu_max,start_index,u_max_eq_index_k);

    % DU constr
    [start_index,du_min_eq_index_k] = expand_equal_index(ndu_min,start_index,du_min_eq_index_k);
    [start_index,du_max_eq_index_k] = expand_equal_index(ndu_max,start_index,du_max_eq_index_k);

    % Y constr
    [start_index,y_min_eq_index_k] = expand_equal_index(ny_min,start_index,y_min_eq_index_k);
    [start_index,y_max_eq_index_k] = expand_equal_index(ny_max,start_index,y_max_eq_index_k);

    % H constr
    [start_index,h_min_eq_index_k] = expand_equal_index(nh_min,start_index,h_min_eq_index_k);
    [start_index,h_max_eq_index_k] = expand_equal_index(nh_max,start_index,h_max_eq_index_k);
end

% S constr
[start_index,s_min_eq_index_k] = expand_equal_index(ns_min,start_index,s_min_eq_index_k);
[start_index,s_max_eq_index_k] = expand_equal_index(ns_max,start_index,s_max_eq_index_k);

mpc.dyn_k = dyn_k';
mpc.su_dyn_k = su_dyn_k';

if mpc.has_s_cnstr
    if mpc.s_cnstr.min_limit
    mpc.s_cnstr.min_eq_index_k = [zeros(mpc.nx,1) s_min_eq_index_k'];
    end
    if mpc.s_cnstr.max_limit
    mpc.s_cnstr.max_eq_index_k = [zeros(mpc.nx,1) s_max_eq_index_k'];
    end
end
if mpc.has_u_cnstr
    if mpc.u_cnstr.min_limit
    mpc.u_cnstr.min_eq_index_k = [u_min_eq_index_k' zeros(mpc.nu,1)];
    end
    if mpc.u_cnstr.max_limit
    mpc.u_cnstr.max_eq_index_k = [u_max_eq_index_k' zeros(mpc.nu,1)];
    end
end
if mpc.has_du_cnstr
    if mpc.du_cnstr.min_limit
    mpc.du_cnstr.min_eq_index_k = [du_min_eq_index_k' zeros(mpc.nu,1)];
    end
    if mpc.du_cnstr.max_limit
    mpc.du_cnstr.max_eq_index_k = [du_max_eq_index_k' zeros(mpc.nu,1)];
    end
end
if mpc.has_y_cnstr
    if mpc.y_cnstr.min_limit
    mpc.y_cnstr.min_eq_index_k = [zeros(mpc.ny,1) y_min_eq_index_k' zeros(mpc.ny,1)];
    end
    if mpc.y_cnstr.max_limit
    mpc.y_cnstr.max_eq_index_k = [zeros(mpc.ny,1) y_max_eq_index_k' zeros(mpc.ny,1)];
    end
end
if mpc.has_h_cnstr
    if mpc.h_cnstr.min_limit
    mpc.h_cnstr.min_eq_index_k = [zeros(mpc.nh-nh_min_0,1) h_min_eq_index_k' zeros(mpc.nh,1)];
    end
    if mpc.h_cnstr.max_limit
    mpc.h_cnstr.max_eq_index_k = [zeros(mpc.nh-nh_max_0,1) h_max_eq_index_k' zeros(mpc.nh,1)];
    end
end

end
