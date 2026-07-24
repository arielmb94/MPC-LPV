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
function mpc = update_mpc_slack_cost(mpc,cnstr,qv_min,qv_max)

if ~isempty(qv_min)
    if isscalar(qv_min) 
        if cnstr.use_k0, mpc.grad_qv_0(cnstr.min_row_v_0) = qv_min; end
        for k = 1:mpc.N-1
            mpc.grad_qv_k(cnstr.min_row_v_k,k) = qv_min;
        end
        if cnstr.use_ter, mpc.grad_qv_ter(cnstr.min_ineqRow_ter) = qv_min; end
    else
        if cnstr.use_k0, mpc.grad_qv_0(cnstr.min_row_v_0) = qv_min(cnstr.rows_k0); end
        for k = 1:mpc.N-1
            mpc.grad_qv_k(cnstr.min_row_v_k,k) = qv_min;
        end
        if cnstr.use_ter, mpc.grad_qv_ter(cnstr.min_ineqRow_ter) = qv_min(cnstr.rows_ter); end
    end
end

if ~isempty(qv_max)
    if isscalar(qv_max) 
        if cnstr.use_k0, mpc.grad_qv_0(cnstr.max_row_v_0) = qv_max; end
        for k = 1:mpc.N-1
            mpc.grad_qv_k(cnstr.max_row_v_k,k) = qv_max;
        end
        if cnstr.use_ter, mpc.grad_qv_ter(cnstr.max_ineqRow_ter) = qv_max; end
    else
        if cnstr.use_k0, mpc.grad_qv_0(cnstr.max_row_v_0) = qv_max(cnstr.rows_k0); end
        for k = 1:mpc.N-1
            mpc.grad_qv_k(cnstr.max_row_v_k,k) = qv_max;
        end
        if cnstr.use_ter, mpc.grad_qv_ter(cnstr.max_ineqRow_ter) = qv_max(cnstr.rows_ter); end
    end
end

end