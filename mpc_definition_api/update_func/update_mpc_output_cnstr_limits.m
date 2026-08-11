% UPDATE_MPC_OUTPUT_CNSTR_LIMITS Update output bounds.
%
%   mpc = UPDATE_MPC_OUTPUT_CNSTR_LIMITS(mpc, y_min, y_max) updates
%
%       y_min_k <= y_k <= y_max_k.
%
%   Use [] to leave either bound unchanged. The corresponding lower or
%   upper bound must first be enabled with INIT_MPC_OUTPUT_CNSTR. Output
%   constraints remain soft after their limits are updated.
%
%   A bound may be a scalar, an ny-by-1 vector, or a time-varying ny-by-L
%   matrix, where L is the number of supplied horizon stages. If L < N, the
%   last supplied column is reused for the remaining stages.
%
%   Inputs:
%     mpc     - Built CHRONOS MPC structure.
%     y_min   - Optional updated lower bound: scalar, ny-by-1, or ny-by-L.
%     y_max   - Optional updated upper bound: scalar, ny-by-1, or ny-by-L.
%
%   Output:
%     mpc     - Updated CHRONOS MPC structure.
%
%   Example - update only the upper output bound:
%
%       mpc = update_mpc_output_cnstr_limits(mpc, [], y_max);
function mpc = update_mpc_output_cnstr_limits(mpc, min, max)

if ~isempty(min)

    if isscalar(min)
        mpc.y_cnstr.min(:, :) = min;

        if mpc.y_cnstr.use_k0
            mpc.y_cnstr.min_0(:) = min;
        end
        if mpc.y_cnstr.use_ter
            mpc.y_cnstr.min_ter(:) = min;
        end
    else

        if mpc.y_cnstr.use_k0
            mpc.y_cnstr.min_0(:) = min(mpc.y_cnstr.rows_k0, 1);
        end

        if size(min, 2) < mpc.N
            mpc.y_cnstr.min(:, :) = fill_vec(mpc.y_cnstr.min, min, 1);
        else
            mpc.y_cnstr.min(:, :) = min(:, 1:mpc.N-1);
        end
        if mpc.y_cnstr.use_ter
            if size(min, 2) < mpc.N
                mpc.y_cnstr.min_ter(:) = min(mpc.y_cnstr.rows_ter, size(min, 2));
            else
                mpc.y_cnstr.min_ter(:) = min(mpc.y_cnstr.rows_ter, mpc.N);
            end
        end
    end
end

if ~isempty(max)

    if isscalar(max)
        mpc.y_cnstr.max(:, :) = max;

        if mpc.y_cnstr.use_k0
            mpc.y_cnstr.max_0(:) = max;
        end
        if mpc.y_cnstr.use_ter
            mpc.y_cnstr.max_ter(:) = max;
        end
    else

        if mpc.y_cnstr.use_k0
            mpc.y_cnstr.max_0(:) = max(mpc.y_cnstr.rows_k0, 1);
        end

        if size(max, 2) < mpc.N
            mpc.y_cnstr.max(:, :) = fill_vec(mpc.y_cnstr.max, max, 1);
        else
            mpc.y_cnstr.max(:, :) = max(:, 1:mpc.N-1);
        end
        if mpc.y_cnstr.use_ter
            if size(max, 2) < mpc.N
                mpc.y_cnstr.max_ter(:) = max(mpc.y_cnstr.rows_ter, size(max, 2));
            else
                mpc.y_cnstr.max_ter(:) = max(mpc.y_cnstr.rows_ter, mpc.N);
            end
        end
    end
end
end
