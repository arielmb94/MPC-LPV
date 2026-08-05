% UPDATE_MIN_SLACK_COST Updates the penalty weights for slack values on the 
% minumum limits of box constraints with soft constraints defined during
% initialization.
%
%   mpc = UPDATE_MIN_SLACK_COST(mpc, cnstr, qv_min) sets to qv_min the
%   penalty weights of the slack values active on the minimum limits of the
%   given constraint and adapts the CHRONOS mpc gradients accordignly. 
%
%   The list of available CHRONOS box constraint structures with slack
%   variables is:
%
%   mpc.s_cnstr            - box constraint on the states
%   mpc.y_cnstr            - box constraint on the output tracking signals
%   mpc.h_cnstr            - box constraint on the user defined constraints
% 
%   INPUTS:
%       mpc                - CHRONOS MPC structure
%       cnstr              - CHRONOS structure for the box constraint to be
%                          updated.                         
%       qv_min             - [ni x 1] or scalar. Penalty weight for 
%                          violating the cnstr minimum soft limits. Higher 
%                          values mean stricter enforcement. ni is the full
%                          size of the constraint (e.g. nx for constraints
%                          on the state vector, nu for constraints on the
%                          control actions and control action rate, ny for
%                          constraints of the tracking signal and nh for
%                          user defined contraints)
%
%   OUTPUTS:
%       mpc                - Updated MPC structure. All necessary 
%                          background math are automatically assembled and 
%                          added to the object.
%
%  USAGE TIPS:
%       - To be used only on constraints that already have soft constraints
%       defined during initialization of the CRHONOS mpc problem
%       - Passing a scalar to qv_min will automatically apply that setting 
%       across all slacks active on the soft constraint.
%       - Passing an [ni x 1] vector you can modify individually the penalty
%       weight for each constraint element.
function mpc = update_mpc_slack_cost(mpc, cnstr, qv_min, qv_max)
% Note: The struct fields cnstr.qv_min and cnstr.qv_max are intentionally 
% NOT updated here to maintain a clean single-output API (returning only mpc). 
% The underlying solver only requires the updated global mpc.grad_qv vectors.

if ~isempty(qv_min)
    if isscalar(qv_min)

        if cnstr.use_k0
            mpc.grad_qv_0(cnstr.min_row_v_0) = qv_min; 
        end
        
        mpc.grad_qv_k(cnstr.min_row_v_k, :) = qv_min;
        
        if cnstr.use_ter
            mpc.grad_qv_ter(cnstr.min_ineqRow_ter) = qv_min; 
        end
    else
        % Initial stage (k = 0)
        if cnstr.use_k0
            mpc.grad_qv_0(cnstr.min_row_v_0) = qv_min(cnstr.rows_k0, 1); 
        end
        
        % Intermediate stages (k = 1 ... N-1)
        if size(qv_min, 2) < mpc.N
            mpc.grad_qv_k(cnstr.min_row_v_k, :) = fill_vec(mpc.grad_qv_k(cnstr.min_row_v_k, :), qv_min, 1);
        else
            mpc.grad_qv_k(cnstr.min_row_v_k, :) = qv_min(:, 1:mpc.N-1);
        end
        
        % Terminal stage (k = N)
        if cnstr.use_ter
            if size(qv_min, 2) < mpc.N
                mpc.grad_qv_ter(cnstr.min_ineqRow_ter) = qv_min(cnstr.rows_ter, end);
            else
                mpc.grad_qv_ter(cnstr.min_ineqRow_ter) = qv_min(cnstr.rows_ter, mpc.N);
            end
        end
    end
end

if ~isempty(qv_max)
    if isscalar(qv_max)
        % Fully vectorized scalar assignment
        if cnstr.use_k0
            mpc.grad_qv_0(cnstr.max_row_v_0) = qv_max; 
        end
        
        mpc.grad_qv_k(cnstr.max_row_v_k, :) = qv_max;
        
        if cnstr.use_ter
            mpc.grad_qv_ter(cnstr.max_ineqRow_ter) = qv_max; 
        end
    else
        % Initial stage (k = 0)
        if cnstr.use_k0
            mpc.grad_qv_0(cnstr.max_row_v_0) = qv_max(cnstr.rows_k0, 1); 
        end
        
        % Intermediate stages (k = 1 ... N-1)
        if size(qv_max, 2) < mpc.N
            mpc.grad_qv_k(cnstr.max_row_v_k, :) = fill_vec(mpc.grad_qv_k(cnstr.max_row_v_k, :), qv_max, 1);
        else
            mpc.grad_qv_k(cnstr.max_row_v_k, :) = qv_max(:, 1:mpc.N-1);
        end
        
        % Terminal stage (k = N)
        if cnstr.use_ter
            if size(qv_max, 2) < mpc.N
                mpc.grad_qv_ter(cnstr.max_ineqRow_ter) = qv_max(cnstr.rows_ter, end);
            else
                mpc.grad_qv_ter(cnstr.max_ineqRow_ter) = qv_max(cnstr.rows_ter, mpc.N);
            end
        end
    end
end
end