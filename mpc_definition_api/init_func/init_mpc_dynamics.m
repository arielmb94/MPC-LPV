%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   mpc = init_mpc_system(mpc,A,B,Bd,C,D,Dd)
%
% Initializes the MPC Discrete-Time dynamical model :
%
%   x+ = A * x + B * u + Bd * d
%
% and the measurement model for the MPC tracking signal:
%
%   y = C * x + D * u + Dd * d
%
% x and u are the state and input vectors, d corresponds to a measured or
% estimated disturbance vector, to be introduced on the appropiate field on
% mpc_solve() during runtime MPC execution.
%
% In:
%   - mpc: CHRONOS mpc structure
%   - A: nx x nx matrix, system matrix
%   - B: nx x nu matrix, input matrix
%   - Bd: nx x nd matrix, disturbance input matrix. If it does not exists,
%   must be set to 0
%   - C: ny x nx matrix, system output matrix
%   - D: ny x nu matrix, output feedtrhough matrix. If it does not exists,
%   must be set to 0
%   - Dd: ny x nd matrix, disturbance output feedtrhough matrix. If it does
%   not exists must be set to 0
%
% Out:
%   - mpc: updated CHRONOS mpc structure
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function mpc = init_mpc_dynamics(mpc,A,B,Bd)
arguments
    mpc
    A = []
    B = []
    Bd = []
end

mpc.Bd = Bd;

%number of states
mpc.nx = size(A,1);  
mpc.A = zeros(mpc.nx,mpc.nx,mpc.N);
mpc.A = fill_mat(mpc.A, A, 1);
%number of control inputs
mpc.nu = size(B,2);  
mpc.B = zeros(mpc.nx,mpc.nu,mpc.N);
mpc.B = fill_mat(mpc.B, B, 1);
%number of disturbance inputs
if any(Bd), mpc.nd = max([size(Bd,2) mpc.nd]); end  
if ~isempty(Bd) && any(Bd(:))
    mpc.dyn_use_d = 1;
    mpc.Bd = zeros(mpc.nx,mpc.nd,mpc.N);
    mpc.Bd = fill_mat(mpc.Bd, Bd, 1);
end

mpc.Nx = mpc.N*mpc.nx;
mpc.Nu = mpc.N*mpc.nu;
mpc.Nd = mpc.N*mpc.nd;

mpc.s = zeros(mpc.nx,mpc.N);
mpc.s_ter = zeros(mpc.nx,1);
mpc.su = zeros(mpc.nu,mpc.N);
mpc.u = zeros(mpc.nu,mpc.N);
mpc.du = zeros(mpc.nu,mpc.N);

if mpc.nd
    mpc.d = zeros(mpc.nd,mpc.N);
end

% Assume C = I*x
C = eye(mpc.nx);
mpc.C = zeros(mpc.nx,mpc.nx,mpc.N-1);
mpc.C = fill_mat(mpc.C, C, 1);
mpc.C_ter = eye(mpc.nx);

mpc.ny = mpc.nx;
mpc.ny_0 = 0;
mpc.ny_ter = mpc.nx;

mpc.y_use_k0 = 0;
mpc.y_rows_k0 = [];
mpc.y_use_ter = 1;
mpc.y_rows_ter = 1:mpc.nx; 

mpc.y_use_s = 1;
mpc.y_use_u = 0;
mpc.y_use_d = 0;

% init y, reference and error vectors
mpc.r_0 = zeros(mpc.ny_0,1);
mpc.y_0 = zeros(mpc.ny_0,1);
mpc.err_0 = zeros(mpc.ny_0,1);
mpc.r = zeros(mpc.ny,mpc.N-1);
mpc.y = zeros(mpc.ny,mpc.N-1);
mpc.err = zeros(mpc.ny,mpc.N-1);
mpc.r_ter = zeros(mpc.ny_ter,1);
mpc.y_ter = zeros(mpc.ny_ter,1);
mpc.err_ter = zeros(mpc.ny_ter,1);

end