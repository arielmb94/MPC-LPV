%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   mpc = update_mpc_Tracking_cost(mpc,Qe)
%
% Updates the tracking error cost weight Qe.
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
function mpc = update_mpc_Tracking_cost(mpc,Qe)

mpc.update_tracking = 1;

len_Q = size(Qe,3);
if len_Q < mpc.N
    mpc.Qe(:,:,:) = fill_mat(mpc.Qe, Qe, 1);
    if mpc.y_use_ter, mpc.Qe_ter(:,:) = mpc.Qe(mpc.y_rows_ter,mpc.y_rows_ter,mpc.N-1); end
else
    mpc.Qe(:,:,:) = Qe(:,:,1:mpc.N-1);
    if mpc.y_use_ter, mpc.Qe_ter(:,:) = Qe(mpc.y_rows_ter,mpc.y_rows_ter,mpc.N); end
end

if mpc.y_use_k0, mpc.Qe_0(:,:) = mpc.Qe(mpc.y_rows_k0,mpc.y_rows_k0,1); end

end