% Initialize random system

sys = drss(5,2,5)
%%
N = 3;              %prediction horizon

A = sys.A;
B = sys.B;
C = sys.C;
D = sys.D;

nx = size(A,1);  %number of states
nd = 1;          %number of distrubance inputs
nu = size(B,2)-1;  %number of control inputs
ny = size(C,1);  %number of measurements

Bd = B(:,nu+nd:end);
B = B(:,1:nu);

Dd = D(:,nu+nd:end);
D = D(:,1:nu);

Nx = (N+1)*nx;
Nu = (N+1)*nu;
Ny = N*ny;
Nd = (N+1)*nd;

x0 = zeros(Nu+Nx,1);
x_prev = rand(nx,1);
u_prev = rand(nu,1);


d = rand(Nd,1);

Qe = diag(100*ones(ny,1))
R = diag(5*ones(nu,1))
P = rand(nx);
P = P*P';
r = ones(ny,1);



% mpc structure

x_min = [];%-10*ones(nx,1);
x_max = [];%10*ones(nx,1);
x_ter_min = -10*ones(nx,1);
x_ter_max = 10*ones(nx,1);
u_min = -10*ones(nu,1);
u_max = 10*ones(nu,1);
du_min = -1*ones(nu,1);
du_max = 1*ones(nu,1);
y_min = -5*ones(ny,1);
y_max = 5*ones(ny,1);


%
mpc = defLtiMpc(N,A,B,C,D,Bd,Dd,Qe,R,x_min,x_max,x_ter_min,x_ter_max,u_min,u_max,du_min,du_max,y_min,y_max)

%%

sigma = 1e-3;
eps_ipopt = 1e-0;

% Compute inequality functions at x0 to compute duality measure
mpc.t = init_t(x0,u_prev,mpc.C,mpc.D,mpc.Dd,d,sigma,mpc.x_min,mpc.x_max,...
    mpc.u_min,mpc.u_max,mpc.du_min,mpc.du_max,mpc.y_min,mpc.y_max, ...
    mpc.N,mpc.Nx,mpc.Nu,mpc.Ny,mpc.nx,mpc.nu,mpc.ny,mpc.nd);

mpc.eta_fwd = 3;
mpc.eta_bck= 4;
mpc.t_max = mpc.m/eps_ipopt;

mpc.max_iter = 1;
