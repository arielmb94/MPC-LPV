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

%% Optimization variebles x
mpc.n = mpc.Nx+(ns_min*mpc.N+ns_max*mpc.N)*2+...
        mpc.Nu+nu_min*mpc.N+nu_max*mpc.N;

x = zeros(mpc.n,1);

u_index = [];
s_index = [];
g_index = [];
v_index = [];

s_index_k = [];
gs_min_index_k = [];
gs_max_index_k = [];
vs_min_index_k = [];
vs_max_index_k = [];
u_index_k = [];
gu_min_index_k = [];
gu_max_index_k = [];

start_index = 1;

% k = 0
dim = [mpc.nu ...      % u
       nu_min nu_max]; % gu
% u
[start_index,u_index,u_index_k] = expand_index(dim,start_index,1,u_index,u_index_k);
% gu
[start_index,g_index,gu_min_index_k] = expand_index(dim,start_index,2,g_index,gu_min_index_k);
[start_index,g_index,gu_max_index_k] = expand_index(dim,start_index,3,g_index,gu_max_index_k);

for k = 1:mpc.N-1

dim = [mpc.nx mpc.nu ...               % s u
       ns_min ns_max nu_min nu_max ... % gs gu
       ns_min ns_max];                 % vs

% s
[start_index,s_index,s_index_k] = expand_index(dim,start_index,1,s_index,s_index_k);

% u
[start_index,u_index,u_index_k] = expand_index(dim,start_index,2,u_index,u_index_k);

% gs
[start_index,g_index,gs_min_index_k] = expand_index(dim,start_index,3,g_index,gs_min_index_k);
[start_index,g_index,gs_max_index_k] = expand_index(dim,start_index,4,g_index,gs_max_index_k);

% gu
[start_index,g_index,gu_min_index_k] = expand_index(dim,start_index,5,g_index,gu_min_index_k);
[start_index,g_index,gu_max_index_k] = expand_index(dim,start_index,6,g_index,gu_max_index_k);

% vs
[start_index,v_index,vs_min_index_k] = expand_index(dim,start_index,7,v_index,vs_min_index_k);
[start_index,v_index,vs_max_index_k] = expand_index(dim,start_index,8,v_index,vs_max_index_k);

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
mpc.u_index = u_index';
mpc.u_index_k = [u_index_k' 0];
mpc.g_index = g_index';
mpc.v_index = v_index';

x(mpc.g_index) = 1;
x(mpc.v_index) = 1;
mpc.slack_index = find(x==1);

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

%% Equality rows 
start_index = 1;

dyn_k = [];
u_min_eq_index_k = [];
u_max_eq_index_k = [];
s_min_eq_index_k = [];
s_max_eq_index_k = [];

% k = 0
% Dynamics
[start_index,dyn_k] = expand_equal_index(mpc.nx,start_index,dyn_k);

% U constr
[start_index,u_min_eq_index_k] = expand_equal_index(nu_min,start_index,u_min_eq_index_k);
[start_index,u_max_eq_index_k] = expand_equal_index(nu_max,start_index,u_max_eq_index_k);

for k = 1:mpc.N-1

    % Dynamics
    [start_index,dyn_k] = expand_equal_index(mpc.nx,start_index,dyn_k);

    % S constr
    [start_index,s_min_eq_index_k] = expand_equal_index(ns_min,start_index,s_min_eq_index_k);
    [start_index,s_max_eq_index_k] = expand_equal_index(ns_max,start_index,s_max_eq_index_k);

    % U constr
    [start_index,u_min_eq_index_k] = expand_equal_index(nu_min,start_index,u_min_eq_index_k);
    [start_index,u_max_eq_index_k] = expand_equal_index(nu_max,start_index,u_max_eq_index_k);
end

% S constr
[start_index,s_min_eq_index_k] = expand_equal_index(ns_min,start_index,s_min_eq_index_k);
[start_index,s_max_eq_index_k] = expand_equal_index(ns_max,start_index,s_max_eq_index_k);

mpc.dyn_k = dyn_k';

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

end
