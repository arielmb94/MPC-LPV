%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   mpc = init_mpc_Tracking_cost(mpc,Qe)
%
% Adds quadratic penalties on the tracking error:
%
%   J += (r - y)' * Qe * (r - y)
%
% y is the tracking feedback signal, defined during the call to 
% init_mpc_system(), the reference vector r is introduced during MPC 
% runtime iterations on the call to mpc_solve().
%
% In:
%   - mpc: CHRONOS mpc structure
%   - Qe: ny x ny square matrix, weights for the quadratic penalty on the
%   tracking error
%
% Out:
%   - mpc: updated CHRONOS mpc structure
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function mpc = init_mpc_Tracking_cost(mpc,Qe)

mpc.tracking_cost = 1;

mpc.Qe = zeros(mpc.ny,mpc.ny,mpc.N-1);
len_Qe = size(Qe,3);
if len_Qe < mpc.N
    mpc.Qe = fill_mat(mpc.Qe, Qe, 1);
    Qe_ter = mpc.Qe(:,:,mpc.N-1);
else
    mpc.Qe = Qe(:,:,1:mpc.N-1);
    Qe_ter = Qe(:,:,mpc.N);
end

if mpc.y_use_k0, mpc.Qe_0 = Qe(mpc.y_rows_k0,mpc.y_rows_k0,1); end
if mpc.y_use_ter, mpc.Qe_ter = Qe_ter(mpc.y_rows_ter,mpc.y_rows_ter); end

end