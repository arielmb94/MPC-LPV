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

%number of measurements
mpc.ny = max([size(C,1) size(D,1)]);  

if ~isempty(C) && any(C(:))
    mpc.y_use_s = 1;

    len_C = size(C,3);
    if len_C < mpc.N
        mpc.C = zeros(mpc.ny,mpc.nx,mpc.N-1);
        mpc.C = fill_mat(mpc.C, C, 1);
        C_ter = mpc.C(:,:,mpc.N-1);
    else
        mpc.C = C(:,:,1:mpc.N-1);
        C_ter = C(:,:,mpc.N);
    end
end

if ~isempty(D) && any(D(:))
    mpc.y_use_u = 1;

    len_D = size(D,3);
    if len_D < mpc.N
        mpc.D = zeros(mpc.ny,mpc.nu,mpc.N-1);
        mpc.D = fill_mat(mpc.D, D, 1);
        D_ter = mpc.D(:,:,mpc.N-1);
    else
        mpc.D = D(:,:,1:mpc.N-1);
        D_ter = D(:,:,mpc.N);
    end
end

if ~isempty(Dd) && any(Dd(:))
    mpc.y_use_d = 1;
    mpc.nd = max([size(Dd,2) mpc.nd]);
    mpc.d = zeros(mpc.nd,mpc.N);

    len_Dd = size(Dd,3);
    if len_Dd < mpc.N
        mpc.Dd = zeros(mpc.ny,mpc.nd,mpc.N-1);
        mpc.Dd = fill_mat(mpc.Dd, Dd, 1);
        Dd_ter = mpc.Dd(:,:,mpc.N-1);
    else
        mpc.Dd = Dd(:,:,1:mpc.N-1);
        Dd_ter = Dd(:,:,mpc.N);
    end
end

% at k = 0, only rows with D!=0 (with dependence on control action u) are
% considered
y_row_0 = find(~all(mpc.D(:,:,1)==0,2));
mpc.ny_0 = length(y_row_0);

if mpc.ny_0
    mpc.y_rows_k0 = y_row_0;
    mpc.y_use_k0 = 1;

    if mpc.y_use_s, mpc.C_0 = mpc.C(mpc.y_rows_k0,:,1); end
    if mpc.y_use_u, mpc.D_0 = mpc.D(mpc.y_rows_k0,:,1); end
    if mpc.y_use_d, mpc.Dd_0 = mpc.Dd(mpc.y_rows_k0,:,1); end
else
    mpc.y_use_k0 = 0;
end

% at k = N, only rows strictly dependent on s are considered
if mpc.y_use_s
    strict_s_rows = any(C_ter~=0,2);
    if mpc.y_use_u, strict_s_rows = strict_s_rows & all(D_ter==0,2); end
    if mpc.y_use_d, strict_s_rows = strict_s_rows & all(Dd_ter==0,2); end

    y_row_ter = find(strict_s_rows);
    mpc.ny_ter = length(y_row_ter);
else
    mpc.ny_ter = 0;
end
if mpc.ny_ter
    mpc.y_rows_ter = y_row_ter;
    mpc.y_use_ter = 1;

    mpc.C_ter = C_ter(mpc.y_rows_ter,:);
else
    mpc.y_rows_ter = [];
    mpc.y_use_ter = 0;

    mpc.C_ter = [];
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