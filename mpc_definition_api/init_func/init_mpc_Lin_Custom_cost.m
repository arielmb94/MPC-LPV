%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   mpc = init_mpc_Lin_Custom_cost(mpc,Cz,Dz,Ddz,Qz,qz)
%
% Adds quadratic and linear penalties on a custom user defined signal z:
%
%   J += z'*Qz*z + qz*z
%
% The signal z is defined as:
%
%   z = Cz * x + Dz * u + Ddz * dz
%
% The user must define the signal z by selecting appropiate values for the
% matrices Cz, Dz, Ddz.
%
% In:
%   - mpc: CHRONOS mpc structure
%   - Cz: nz x nx matrix, states output matrix
%   - Dz: nz x nu matrix, input feedtrhough matrix
%   - Ddz: nz x ndi matrix, disturbance feedtrhough matrix
%   - Qz (optional): nz x nz square matrix, weights for the quadratic
%   penalty term on the user defined signal z
%   - qz (optional): nz column vector, weights for the linear penalty term
%   on the user defined signal z. Use the linear penalty term only in the
%   case that z takes only positive values
%
% Out:
%   - mpc: updated CHRONOS mpc structure
%
% Example:
% We want to penalize the variation of the MPC control action with respect 
% a given value u_star, e.g. we want to minimize z = u - u_star.
% To achieve this we define z by selecting:
%   - Cz = [0 0 ... 0]
%   - Dz = [1]
%   - Ddz = [-1]
% which corresponds to z = [0 0 ... 0] * x + [1] * u + [-1] * di
% The "disturbance" term on z corresponds to u_star, to be introduced on 
% the appropiate field on mpc_solve() during runtime MPC execution.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function mpc = init_mpc_Lin_Custom_cost(mpc,Cz,Dz,Dsuz,Ddz,Qz,qz)
arguments
    mpc
    Cz = []
    Dz = []
    Dsuz = []
    Ddz = []
    Qz = [];
    qz = [];
end

mpc.Qz = Qz;
mpc.qz = qz;

mpc.Cz = Cz;
mpc.Dz = Dz;
mpc.Dsuz = Dsuz;
mpc.Ddz = Ddz;

%number of general inequalities
if ~isempty(Cz) && max(any(Cz))
    mpc.nz = size(Cz,1);  
elseif ~isempty(Dz) && max(any(Dz))
    mpc.nz = size(Dz,1);
elseif ~isempty(Dsuz) && max(any(Dsuz))
    mpc.nz = size(Dsuz,1);
end
mpc.ndz = size(Ddz,2);  %number of disturbance inputs to general inequalities

if ~isempty(Dsuz) && max(any(Dsuz)), mpc.has_du = 1; end

if ~isempty(mpc.Cz) && max(any(mpc.Cz))
    mpc.z_use_s = 1;
end
if ~isempty(mpc.Dz) && max(any(mpc.Dz))
    mpc.z_use_u = 1;
end
if ~isempty(mpc.Dsuz) && max(any(mpc.Dsuz))
    mpc.z_use_su = 1;
end
if ~isempty(mpc.Ddz) && max(any(mpc.Ddz))
    mpc.z_use_d = 1;
end

mpc.z_use_k0 = 0;
mpc.z_use_ter = 0;
% at k = 0, only rows with Dz!=0 (with dependence on control action u) are
% considered
z_row_0 = find(~all(Dz==0,2));
mpc.nz_0 = length(z_row_0);

if mpc.nz_0
    mpc.z_rows_k0 = z_row_0;
    mpc.z_use_k0 = 1;

    if mpc.z_use_s, mpc.Cz_0 = Cz(mpc.z_rows_k0,:); end
    if mpc.z_use_u, mpc.Dz_0 = Dz(mpc.z_rows_k0,:); end
    if mpc.z_use_su, mpc.Dsuz_0 = Dsuz(mpc.z_rows_k0,:); end
    if mpc.z_use_d, mpc.Ddz_0 = Ddz(mpc.z_rows_k0,:); end
end

% at k = N, only rows strictly dependent on s are considered
if  ~isempty(Cz)
    strict_s_rows = any(Cz~=0,2);
    if mpc.z_use_u, strict_s_rows = strict_s_rows & all(Dz==0,2); end
    if mpc.z_use_su, strict_s_rows = strict_s_rows & all(Dsuz==0,2); end
    if mpc.z_use_d, strict_s_rows = strict_s_rows & all(Ddz==0,2); end

    z_row_ter = find(strict_s_rows);
    mpc.nz_ter = length(z_row_ter);
else
    mpc.nz_ter = 0;
end
if mpc.nz_ter
    mpc.z_rows_ter = z_row_ter;
    mpc.z_use_ter = 1;

    mpc.Cz_ter = Cz(mpc.z_rows_ter,:); 
end

% init z vector
if mpc.z_use_k0, mpc.z_0 = zeros(mpc.nz_0,1); else, mpc.z_0=[]; end
mpc.z = zeros(mpc.nz,mpc.N-1);
if mpc.z_use_ter, mpc.z_ter = zeros(mpc.nz_ter,1); else, mpc.z_ter=[]; end
% init d vector
if mpc.ndz
    mpc.dz = zeros(mpc.ndz,mpc.N);
end

if ~isempty(Qz)
    mpc.quad_custom_cost = 1;

    if mpc.z_use_k0, mpc.Qz_0 = Qz(mpc.z_rows_k0,mpc.z_rows_k0); end
    if mpc.z_use_ter, mpc.Qz_ter = Qz(mpc.z_rows_ter,mpc.z_rows_ter); end
end
if ~isempty(qz)
    mpc.lin_custom_cost = 1;

    if mpc.z_use_k0, mpc.qz_0 = qz(mpc.z_rows_k0); end
    if mpc.z_use_ter, mpc.qz_ter = qz(mpc.z_rows_ter); end
end
 
end