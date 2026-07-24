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
function mpc = init_mpc_output(mpc,C,D,Dd)
arguments
    mpc
    C = []
    D = []
    Dd = []
end

mpc.C = C;
mpc.D = D;
mpc.Dd = Dd;

mpc.nx = size(mpc.A,1);  %number of states
mpc.nu = size(mpc.B,2);  %number of control inputs
if any(Dd), mpc.nd = max([size(mpc.Dd,2) mpc.nd]); end  %number of disturbance inputs
mpc.ny = max([size(mpc.C,1) size(mpc.D,1)]);  %number of measurements

if mpc.nd
    mpc.d = zeros(mpc.nd,mpc.N);
end

if ~isempty(mpc.C) && max(any(mpc.C))
    mpc.y_use_s = 1;
end
if ~isempty(mpc.D) && max(any(mpc.D))
    mpc.y_use_u = 1;
end
if ~isempty(mpc.Dd) && max(any(mpc.Dd))
    mpc.y_use_d = 1;
end

% at k = 0, only rows with D!=0 (with dependence on control action u) are
% considered
y_row_0 = find(~all(D==0,2));
mpc.ny_0 = length(y_row_0);

if mpc.ny_0
    mpc.y_rows_k0 = y_row_0;
    mpc.y_use_k0 = 1;

    if mpc.y_use_s, mpc.C_0 = C(mpc.y_rows_k0,:); end
    if mpc.y_use_u, mpc.D_0 = D(mpc.y_rows_k0,:); end
    if mpc.y_use_d, mpc.Dd_0 = Dd(mpc.y_rows_k0,:); end
else
    mpc.y_use_k0 = 0;
end

% at k = N, only rows strictly dependent on s are considered
if ~isempty(C)
    strict_s_rows = any(C~=0,2);
    if mpc.y_use_u, strict_s_rows = strict_s_rows & all(D==0,2); end
    if mpc.y_use_d, strict_s_rows = strict_s_rows & all(Dd==0,2); end

    y_row_ter = find(strict_s_rows);
    mpc.ny_ter = length(y_row_ter);
else
    mpc.ny_ter = 0;
end
if mpc.ny_ter
    mpc.y_rows_ter = y_row_ter;
    mpc.y_use_ter = 1;

    mpc.C_ter = C(mpc.y_rows_ter,:);
end

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