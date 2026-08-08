% UPDATE_MPC_TRACKING_COST Update the tracking-error weight.
%
%   mpc = UPDATE_MPC_TRACKING_COST(mpc, Qe) replaces the weight in
%
%       J_tracking += 0.5*(r_k-y_k)'*Qe_k*(r_k-y_k).
%
%   Qe may be ny-by-ny or ny-by-ny-by-L, where L is the number of supplied
%   horizon stages. If L < N, the last supplied stage is reused for the
%   remaining stages.
%
%   Inputs:
%     mpc     - Built CHRONOS MPC structure.
%     Qe      - Updated tracking weight, size ny-by-ny or ny-by-ny-by-L.
%
%   Output:
%     mpc     - Updated CHRONOS MPC structure.
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
