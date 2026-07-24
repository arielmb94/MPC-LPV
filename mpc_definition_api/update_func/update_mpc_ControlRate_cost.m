%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   mpc = update_mpc_DiffControl_cost(mpc,Rdu)
%
% Modifies the weight Rdu for the quadratic penalty term on the control 
% action variation between sampling instances. The function then updates 
% the MPC gradients and Hessians accordingly.
%
% In:
%   - mpc: CHRONOS mpc structure
%   - Rdu: nu x nu square matrix, weights for the quadratic penalty term on
%   the control control variation between sampling instances
%
% Out:
%   - mpc: updated CHRONOS mpc structure
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function mpc = update_mpc_ControlRate_cost(mpc,Rdu)

mpc.recompute_cost_hess = 1;
mpc.Rdu(:,:) = Rdu;

for k = 1:mpc.N-1
    mpc.gradRateCtrl_Rdu_k(:,:,k) = [-Rdu;Rdu];
    mpc.H_RateCtrl_k(:,:,k) = [Rdu -Rdu;-Rdu Rdu];
end

end