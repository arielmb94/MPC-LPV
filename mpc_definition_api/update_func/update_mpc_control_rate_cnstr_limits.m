% UPDATE_MPC_CONTROL_RATE_CNSTR_LIMITS Update control-rate bounds.
%
%   mpc = UPDATE_MPC_CONTROL_RATE_CNSTR_LIMITS(mpc, du_min, du_max) updates
%
%       du_min_k <= Delta_u_k <= du_max_k.
%
%   Use [] to leave either bound unchanged. The corresponding lower or
%   upper bound must first be enabled with INIT_MPC_CONTROL_RATE_CNSTR.
%
%   A bound may be a scalar, an nu-by-1 vector, or a time-varying nu-by-L
%   matrix, where L is the number of supplied horizon stages. If L < N, the
%   last supplied column is reused for the remaining stages.
%
%   Inputs:
%     mpc     - Built CHRONOS MPC structure.
%     du_min  - Optional updated lower bound: scalar, nu-by-1, or nu-by-L.
%     du_max  - Optional updated upper bound: scalar, nu-by-1, or nu-by-L.
%
%   Output:
%     mpc     - Updated CHRONOS MPC structure.
%
%   Example - update only the rate lower bound:
%
%       mpc = update_mpc_control_rate_cnstr_limits(mpc, du_min, []);
function mpc = update_mpc_control_rate_cnstr_limits(mpc,min,max)

if ~isempty(min)

    if isscalar(min)
        mpc.du_cnstr.min(:,:) = min;
    elseif size(min,2) < mpc.N
        mpc.du_cnstr.min(:,:) = fill_vec(mpc.du_cnstr.min, min, 1);
    else
        mpc.du_cnstr.min(:,:) = min(:,1:mpc.N);
    end
    
end

if ~isempty(max)

    if isscalar(max)
        mpc.du_cnstr.max(:,:) = max;
    elseif size(max,2) < mpc.N
        mpc.du_cnstr.max(:,:) = fill_vec(mpc.du_cnstr.max, max, 1);
    else
        mpc.du_cnstr.max(:,:) = max(:,1:mpc.N);
    end

end

end
