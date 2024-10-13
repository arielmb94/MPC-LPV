%% Parameters

Ab = 1;
g = 9.81;

h1 = 0.45;
h2 = 0.45;

Ts = 0.01;

%% LTI system

A = [-sqrt(2*g)*sqrt(h1)/(Ab*h1) 0 0;
     sqrt(2*g)*sqrt(h1)/(Ab*h1) -sqrt(2*g)*sqrt(h2)/(Ab*h2) 0;
     0 -1 0];

B = [1/Ab; 0; 0];

Bd = [0 0 1]';

C = [0 1 0;
     0 0 1];

sys_ct = ss(A,[B Bd],C,0);
sys = c2d(sys_ct,Ts)

N = 25;              %prediction horizon

A = sys.A;
B = sys.B(:,1);
Bd = sys.B(:,2);
C = sys.C;
D = sys.D(:,1);
Dd = [];
%%

nx = length(A);  % number of states
nu = size(B,2);  % number of control inputs
ny = size(C,1);  % number of measurements
nd = size(Bd,2);  % number of distrubance inputs

Nx = (N+1)*nx;
Nu = (N+1)*nu;
Ny = N*ny;
Nd = (N+1)*nd;

x0 = 0.45*ones(Nu+Nx,1);
x_prev = [h1; h2; 0];
u_prev = zeros(nu,1);

Qe = [1 0; 0 15];
R = diag(1*ones(nu,1));

% mpc structure

x_min = [0.05 0.05 -100]';
x_max = [1 1 100]';
x_ter_min = [];%[0.05 0.05 -100]';
x_ter_max = [];%[1 1 100]';
u_min = 0*ones(nu,1);
u_max = 10*ones(nu,1);
du_min = -1*ones(nu,1);
du_max = 1*ones(nu,1);
y_min = [];
y_max = [];

%
mpc = defLtiMpc(N,A,B,C,D,Bd,Dd,Qe,R,x_min,x_max,x_ter_min,x_ter_max,u_min,u_max,du_min,du_max,y_min,y_max)

%%

sigma = 1e-3;
eps_ipopt = 2e-0;

% Compute inequality functions at x0 to compute duality measure
mpc.t = init_t(x0,u_prev,mpc.C,mpc.D,mpc.Dd,[],sigma,mpc.x_min,mpc.x_max,...
    mpc.u_min,mpc.u_max,mpc.du_min,mpc.du_max,mpc.y_min,mpc.y_max, ...
    mpc.N,mpc.Nx,mpc.Nu,mpc.Ny,mpc.Nd,mpc.nx,mpc.nu,mpc.ny,mpc.nd);

mpc.eta_fwd = 4;
mpc.eta_bck= 10;
mpc.t_max = mpc.m/eps_ipopt;

mpc.max_iter = 1;


