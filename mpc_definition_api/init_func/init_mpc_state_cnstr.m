% INIT_MPC_STATE_CNSTR Defines state constraints and soft-constraint penalties.
%
%   mpc = INIT_MPC_STATE_CNSTR(mpc, x_min, x_max) sets strict (hard) lower and 
%   upper bounds on the state variables. The solver will strictly enforce these 
%   limits. This is suitable for absolute physical boundaries, but may cause 
%   the solver to crash (go infeasible) if a disturbance pushes the system too far.
%
%   mpc = INIT_MPC_STATE_CNSTR(mpc, x_min, x_max, x_min_slack_active, x_max_slack_active, qv_min, qv_max) 
%   allows you to define specific bounds as "soft" constraints. Soft constraints 
%   can be safely violated during massive disturbances to keep the solver running, 
%   while applying a customizable penalty to drive the state back within limits 
%   as quickly as possible.
%
%   INPUTS:
%       mpc                - CHRONOS MPC structure.
%       x_min              - [nx x 1] Array of lower state limits (use [] if none).
%       x_max              - [nx x 1] Array of upper state limits (use [] if none).
%       qv_min             - (Optional) [nx x 1] or scalar. Penalty weight for violating 
%                            the x_min soft limits. Higher values mean stricter enforcement.
%       qv_max             - (Optional) [nx x 1] or scalar. Penalty weight for violating 
%                            the x_max soft limits. Higher values mean stricter enforcement.
%
%   OUTPUTS:
%       mpc                - Updated MPC structure. All necessary background math 
%                            (constraint gradients, Hessians, and slack variables) 
%                            are automatically assembled and added to the object.
%
%   USAGE TIPS:
%       - If qv_min or qv_max are not passed, the soft constraint penalty 
%         weight will default to the value stored in mpc.qv
%       - Passing a scalar to qv inputs will automatically apply 
%         that setting across all constrained states.
function mpc = init_mpc_state_cnstr(mpc,x_min,x_max,qv_min,qv_max)
arguments
    mpc
    x_min = []
    x_max = []
    qv_min = []
    qv_max = []
end

% INPUT DIMENSION VALIDATION 
validate_column_vector(x_min, mpc.nx, 'x_min');
validate_column_vector(x_max, mpc.nx, 'x_max');
validate_column_vector(qv_min, mpc.nx, 'qv_min');
validate_column_vector(qv_max, mpc.nx, 'qv_max');

% Expand scalars to full local vectors if needed
if isscalar(x_min), x_min = x_min * ones(mpc.nx, 1); end
if isscalar(x_max), x_max = x_max * ones(mpc.nx, 1); end

s_cnstr.min = x_min;
s_cnstr.max = x_max;
s_cnstr.use_k0 = 0;
s_cnstr.use_ter = 0;
s_cnstr.rows_k0 = [];
s_cnstr.rows_ter = [];

mpc.has_s_cnstr = 1;

if ~isempty(s_cnstr.min)

    s_cnstr.min_limit = 1;
   
    mpc.ng_k(2:3) = mpc.ng_k(2:3) + mpc.nx;
    mpc.nv_k(2:3) = mpc.nv_k(2:3) + mpc.nx;

    s_cnstr.g_min_index_k = [];
    s_cnstr.v_min_index_k = [];

    % Initialize Penalty term for new slack variables
    if isempty(qv_min)
        % if qv isnt defined, it is not initialized until build_chronos_mpc(),
        % but we need to make space
        qv_min = zeros(mpc.nx,1);
    elseif length(qv_min) == 1
        qv_min = qv_min*ones(mpc.nx,1);
    end

    s_cnstr.qv_min = qv_min;

else
    s_cnstr.min_limit = 0;
end

if ~isempty(s_cnstr.max)

    s_cnstr.max_limit = 1;

    mpc.ng_k(2:3) = mpc.ng_k(2:3) + mpc.nx;
    mpc.nv_k(2:3) = mpc.nv_k(2:3) + mpc.nx;

    s_cnstr.g_max_index_k = [];
    s_cnstr.v_max_index_k = [];

    % Initialize Penalty term for new slack variables
    if isempty(qv_max)
        % if qv isnt defined, it is not initialized until build_chronos_mpc(),
        % but we need to make space
        qv_max = zeros(mpc.nx,1);
    elseif length(qv_max) == 1
        qv_max = qv_max*ones(mpc.nx,1);
    end

    s_cnstr.qv_max = qv_max;

else
    s_cnstr.max_limit = 0; 
end

mpc.s_cnstr = s_cnstr;

end