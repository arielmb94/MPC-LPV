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
function mpc = update_mpc_output_cnstr_limits(mpc, min, max)

if ~isempty(min)

    if isscalar(min)
        mpc.y_cnstr.min(:, :) = min;
        mpc.bi_k(mpc.y_cnstr.min_ineqRow_k, :) = -min;

        if mpc.y_cnstr.use_k0
            mpc.y_cnstr.min_0(:) = min;
            mpc.bi_0(mpc.y_cnstr.min_ineqRow_0) = -min;
        end
        if mpc.y_cnstr.use_ter
            mpc.y_cnstr.min_ter(:) = min;
            mpc.bi_ter(mpc.y_cnstr.min_ineqRow_ter) = -min;
        end
    else

        if mpc.y_cnstr.use_k0
            mpc.y_cnstr.min_0(:) = min(mpc.y_cnstr.rows_k0, 1);
            mpc.bi_0(mpc.y_cnstr.min_ineqRow_0) = -mpc.y_cnstr.min_0;
        end

        if size(min, 2) < mpc.N
            mpc.y_cnstr.min(:, :) = fill_vec(mpc.y_cnstr.min, min, 1);
        else
            mpc.y_cnstr.min(:, :) = min(:, 1:mpc.N-1);
        end
        mpc.bi_k(mpc.y_cnstr.min_ineqRow_k, :) = -mpc.y_cnstr.min;        

        if mpc.y_cnstr.use_ter
            if size(min, 2) < mpc.N
                mpc.y_cnstr.min_ter(:) = min(mpc.y_cnstr.rows_ter, size(min, 2));
            else
                mpc.y_cnstr.min_ter(:) = min(mpc.y_cnstr.rows_ter, mpc.N);
            end
            mpc.bi_ter(mpc.y_cnstr.min_ineqRow_ter) = -mpc.y_cnstr.min_ter;
        end
    end
end

if ~isempty(max)

    if isscalar(max)
        mpc.y_cnstr.max(:, :) = max;
        mpc.bi_k(mpc.y_cnstr.max_ineqRow_k, :) = max;

        if mpc.y_cnstr.use_k0
            mpc.y_cnstr.max_0(:) = max;
            mpc.bi_0(mpc.y_cnstr.max_ineqRow_0) = max;
        end
        if mpc.y_cnstr.use_ter
            mpc.y_cnstr.max_ter(:) = max;
            mpc.bi_ter(mpc.y_cnstr.max_ineqRow_ter) = max;
        end
    else

        if mpc.y_cnstr.use_k0
            mpc.y_cnstr.max_0(:) = max(mpc.y_cnstr.rows_k0, 1);
            mpc.bi_0(mpc.y_cnstr.max_ineqRow_0) = mpc.y_cnstr.max_0;
        end

        if size(max, 2) < mpc.N
            mpc.y_cnstr.max(:, :) = fill_vec(mpc.y_cnstr.max, max, 1);
        else
            mpc.y_cnstr.max(:, :) = max(:, 1:mpc.N-1);
        end
        mpc.bi_k(mpc.y_cnstr.max_ineqRow_k, :) = mpc.y_cnstr.max;

        if mpc.y_cnstr.use_ter
            if size(max, 2) < mpc.N
                mpc.y_cnstr.max_ter(:) = max(mpc.y_cnstr.rows_ter, size(max, 2));
            else
                mpc.y_cnstr.max_ter(:) = max(mpc.y_cnstr.rows_ter, mpc.N);
            end
            mpc.bi_ter(mpc.y_cnstr.max_ineqRow_ter) = mpc.y_cnstr.max_ter;
        end
    end
end
end