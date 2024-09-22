% Initialize random system

sys = drss(5,1,2)
%%
N = 10;              %prediction horizon

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

x_min = -10*ones(nx,1);
x_max = 10*ones(nx,1);
u_min = -10*ones(nu,1);
u_max = 10*ones(nu,1);
du_min = -1*ones(nu,1);
du_max = 1*ones(nu,1);
y_min = -5*ones(ny,1);
y_max = 5*ones(ny,1);



%
mpc = defLtiMpc(N,A,B,C,D,Bd,Dd,Qe,R,x_min,x_max,u_min,u_max,du_min,du_max,y_min,y_max)
