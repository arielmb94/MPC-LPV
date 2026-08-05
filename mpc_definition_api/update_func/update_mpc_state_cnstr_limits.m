%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   mpc = update_mpc_u_cnstr(mpc,u_min,u_max)
%
% Modifies the constraints limits on the control action
%
% In:
%   - mpc: CHRONOS mpc structure
%   - u_min (optional): nu column vector, lower bound constraint values on
%   the control action
%   - u_max (optional): nu column vector, upper bound constraint values on
%   the control action
%
% Out:
%   - mpc: updated CHRONOS mpc structure
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function mpc = update_mpc_state_cnstr_limits(mpc,min,max)

if ~isempty(min)

    if isscalar(min)
        mpc.s_cnstr.min(:,:) = min;
    elseif size(min,2) < mpc.N
        mpc.s_cnstr.min(:,:) = fill_vec(mpc.s_cnstr.min, min, 1);
    else
        mpc.s_cnstr.min(:,:) = min(:,1:mpc.N);
    end

    mpc.bi_k(mpc.s_cnstr.min_ineqRow_k,:) = -mpc.s_cnstr.min(:,1:mpc.N-1);
    mpc.bi_ter(mpc.s_cnstr.min_ineqRow_ter) = -mpc.s_cnstr.min(:,mpc.N);
end

if ~isempty(max)
    
    if isscalar(max)
        mpc.s_cnstr.max(:,:) = max;
    elseif size(max,2) < mpc.N
        mpc.s_cnstr.max(:,:) = fill_vec(mpc.s_cnstr.max, max, 1);
    else
        mpc.s_cnstr.max(:,:) = max(:,1:mpc.N);
    end

    mpc.bi_k(mpc.s_cnstr.max_ineqRow_k,:) = mpc.s_cnstr.max(:,1:mpc.N-1);
    mpc.bi_ter(mpc.s_cnstr.max_ineqRow_ter) = mpc.s_cnstr.max(:,mpc.N);
end

end