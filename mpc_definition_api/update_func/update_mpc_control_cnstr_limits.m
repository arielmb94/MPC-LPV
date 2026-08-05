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
function mpc = update_mpc_control_cnstr_limits(mpc,min,max)

if ~isempty(min)

    if isscalar(min)
        mpc.u_cnstr.min(:,:) = min;
    elseif size(min,2) < mpc.N
        mpc.u_cnstr.min(:,:) = fill_vec(mpc.u_cnstr.min, min, 1);
    else
        mpc.u_cnstr.min(:,:) = min(:,1:mpc.N);
    end
    
    mpc.bi_0(mpc.u_cnstr.min_ineqRow_0) = -mpc.u_cnstr.min(:,1);
    mpc.bi_k(mpc.u_cnstr.min_ineqRow_k,:) = -mpc.u_cnstr.min(:,2:mpc.N);
end

if ~isempty(max)

    if isscalar(max)
        mpc.u_cnstr.max(:,:) = max;
    elseif size(max,2) < mpc.N
        mpc.u_cnstr.max(:,:) = fill_vec(mpc.u_cnstr.max, max, 1);
    else
        mpc.u_cnstr.max(:,:) = max(:,1:mpc.N);
    end

    mpc.bi_0(mpc.u_cnstr.max_ineqRow_0) = mpc.u_cnstr.max(:,1);
    mpc.bi_k(mpc.u_cnstr.max_ineqRow_k,:) = mpc.u_cnstr.max(:,2:mpc.N);
end

end