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

mpc.Qe(:,:) = Qe;   
if mpc.y_use_k0, mpc.Qe_0(:,:) = Qe(mpc.y_rows_k0,mpc.y_rows_k0); end
if mpc.y_use_ter, mpc.Qe_ter(:,:) = Qe(mpc.y_rows_ter,mpc.y_rows_ter); end

end