%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   mpc = update_mpc_Control_cost(mpc,Ru,ru)
%
% Modifies the weights Ru and ru for the quadratic and linear control 
% penalty terms on the control action. The function then updates the MPC 
% gradients and Hessians accordingly.
%
% In:
%   - mpc: CHRONOS mpc structure
%   - Ru (optional): nu x nu square matrix, weights for the quadratic
%   penalty term on the control action.
%   - ru (optional): nu column vector, weights for the linear penalty term
%   on the control action. IMPORTANT: Use linear penalties only in the case
%   that the control action takes positive values only.
%
% Out:
%   - mpc: updated CHRONOS mpc structure
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function mpc = update_mpc_Control_cost(mpc,Ru,ru)

if ~isempty(Ru)
    mpc.recompute_cost_hess = 1;

    len_R = size(Ru,3);
    if len_R < mpc.N
        mpc.Ru(:,:,:) = fill_mat(mpc.Ru, Ru, 1);
    else
        mpc.Ru(:,:,:) = Ru(:,:,1:mpc.N);
    end
end

if ~isempty(ru)

    len_r = size(ru,3);
    if len_r < mpc.N
        mpc.ru(:,:) = fill_vec(mpc.ru, ru, 1);
    else
        mpc.ru(:,:) = ru(:,:,1:mpc.N);
    end
end

end