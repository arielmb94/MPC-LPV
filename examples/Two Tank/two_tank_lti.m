%% Parameters

Ab = 1;
g = 9.81;

h1 = 0.45;
h2 = 0.45;

Ts = 0.01;

%% LTI system

A = [-sqrt(2*g)*sqrt(h1)/(Ab*h1) 0;
     sqrt(2*g)*sqrt(h1)/(Ab*h1) -sqrt(2*g)*sqrt(h2)/(Ab*h2)];

B = [1/Ab; 0];

C = [0 1];
%C = eye(2);

sys_ct = ss(A,B,C,0);
sys = c2d(sys_ct,Ts)

N = 2;              %prediction horizon

A = sys.A;
B = sys.B;
Bd = [];
C = sys.C;
D = sys.D;
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
x_prev = [h1; h2];
u_prev = 0.45;

Qe = diag(30*ones(ny,1));
R = diag(1*ones(nu,1));

% mpc structure

x_min = 0.05*ones(nx,1);
x_max = 1*ones(nx,1);
x_ter_min = [];%;0.05*ones(nx,1);
x_ter_max = [];%1*ones(nx,1);
u_min = 0*ones(nu,1);
u_max = 10*ones(nu,1);
du_min = -0.1*ones(nu,1);
du_max = 0.1*ones(nu,1);
y_min = [];
y_max = [];

%
mpc = defLtiMpc(N,A,B,C,D,Bd,Dd,Qe,R,x_min,x_max,x_ter_min,x_ter_max,u_min,u_max,du_min,du_max,y_min,y_max)

%%

mpc.t = 50;

mpc.Beta = 0.75;
mpc.min_l = 0.99;

%%
mpc.ter_ingredients = 1;
mpc.ter_constraint = 0;
mpc.x_ref_is_y = 0;

[K,S] = dlqr(A,B,diag([30 30]),R)

mpc.P = S;
mpc.hessTerminalCost = zeros(Nx+Nu,Nx+Nu);
mpc.hessTerminalCost(nx*(N) + nu*(N+1) + 1 : end,nx*(N) + nu*(N+1) + 1 : end) = ...
    S;
