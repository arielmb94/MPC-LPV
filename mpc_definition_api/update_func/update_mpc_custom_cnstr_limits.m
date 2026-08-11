% UPDATE_MPC_CUSTOM_CNSTR_LIMITS Update custom-signal bounds.
%
%   mpc = UPDATE_MPC_CUSTOM_CNSTR_LIMITS(mpc, h_min, h_max) updates
%
%       h_min_k <= h_k <= h_max_k.
%
%   Use [] to leave either bound unchanged. The corresponding lower or
%   upper bound must first be enabled with INIT_MPC_CUSTOM_CNSTR. To
%   update the definition of h_k, use UPDATE_MPC_CUSTOM_CNSTR_VECTOR.
%
%   A bound may be a scalar, an nh-by-1 vector, or a time-varying nh-by-L
%   matrix, where L is the number of supplied horizon stages. If L < N, the
%   last supplied column is reused for the remaining stages.
%
%   Inputs:
%     mpc     - Built CHRONOS MPC structure.
%     h_min   - Optional updated lower bound: scalar, nh-by-1, or nh-by-L.
%     h_max   - Optional updated upper bound: scalar, nh-by-1, or nh-by-L.
%
%   Output:
%     mpc     - Updated CHRONOS MPC structure.
%
%   Example - update both custom bounds:
%
%       mpc = update_mpc_custom_cnstr_limits(mpc, h_min, h_max);
function mpc = update_mpc_custom_cnstr_limits(mpc, min, max)

if ~isempty(min)
    if isscalar(min)
        % In-place scalar expansion across pre-allocated buffers
        mpc.h_cnstr.min(:, :) = min;

        if mpc.h_cnstr.use_k0
            mpc.h_cnstr.min_0(:) = min;
        end
        if mpc.h_cnstr.use_ter
            mpc.h_cnstr.min_ter(:) = min;
        end
    else

        if mpc.h_cnstr.use_k0
            mpc.h_cnstr.min_0(:) = min(mpc.h_cnstr.rows_k0, 1);
        end

        if size(min, 2) < mpc.N
            mpc.h_cnstr.min(:, :) = fill_vec(mpc.h_cnstr.min, min, 1);
        else
            mpc.h_cnstr.min(:, :) = min(:, 1:mpc.N-1);
        end
        if mpc.h_cnstr.use_ter
            if size(min, 2) < mpc.N
                mpc.h_cnstr.min_ter(:) = min(mpc.h_cnstr.rows_ter, size(min, 2));
            else
                mpc.h_cnstr.min_ter(:) = min(mpc.h_cnstr.rows_ter, mpc.N);
            end
        end
    end
end

if ~isempty(max)
    if isscalar(max)
        mpc.h_cnstr.max(:, :) = max;

        if mpc.h_cnstr.use_k0
            mpc.h_cnstr.max_0(:) = max;
        end
        if mpc.h_cnstr.use_ter
            mpc.h_cnstr.max_ter(:) = max;
        end
    else

        if mpc.h_cnstr.use_k0
            mpc.h_cnstr.max_0(:) = max(mpc.h_cnstr.rows_k0, 1);
        end

        if size(max, 2) < mpc.N
            mpc.h_cnstr.max(:, :) = fill_vec(mpc.h_cnstr.max, max, 1);
        else
            mpc.h_cnstr.max(:, :) = max(:, 1:mpc.N-1);
        end
        if mpc.h_cnstr.use_ter
            if size(max, 2) < mpc.N
                mpc.h_cnstr.max_ter(:) = max(mpc.h_cnstr.rows_ter, size(max, 2));
            else
                mpc.h_cnstr.max_ter(:) = max(mpc.h_cnstr.rows_ter, mpc.N);
            end
        end
    end
end
end
