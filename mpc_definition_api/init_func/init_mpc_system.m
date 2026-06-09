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
function mpc = init_mpc_system(mpc,A,B,Bd,C,D,Dd)

mpc.A = A;
mpc.B = B;
if max(any(Bd))
    mpc.Bd = Bd;
    mpc.dyn_use_d = 1;
end
mpc.C = C;
if max(any(D))
    mpc.D = D;
end
if max(any(Dd))
    mpc.Dd = Dd;
end

mpc.nx = size(mpc.A,1);  %number of states
mpc.nu = size(mpc.B,2);  %number of control inputs
mpc.nd = size(mpc.Bd,2);  %number of disturbance inputs
mpc.ny = size(mpc.C,1);  %number of measurements

mpc.Nx = mpc.N*mpc.nx;
mpc.Nu = mpc.N*mpc.nu;
mpc.Nd = mpc.N*mpc.nd;

mpc.s = zeros(mpc.Nx,1);
mpc.s_all = zeros(mpc.Nx+mpc.nx,1);
mpc.s_ter = zeros(mpc.nx,1);
mpc.u = zeros(mpc.Nu,1);
mpc.du = zeros(mpc.Nu,1);


if mpc.nd
    mpc.d = zeros(mpc.nd,mpc.N);
end

if mpc.ny
    mpc.Ny = (mpc.N-1)*mpc.ny;

    mpc.r = zeros(mpc.ny,mpc.N-1);
    mpc.y = zeros(mpc.ny,mpc.N-1);
    mpc.err = zeros(mpc.ny,mpc.N-1);
end

end