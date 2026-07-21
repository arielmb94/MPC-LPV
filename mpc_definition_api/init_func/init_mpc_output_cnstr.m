% INIT_MPC_OUTPUT_CNSTR Defines output constraints and soft-constraint penalties.
%
%   mpc = INIT_MPC_OUTPUT_CNSTR(mpc, y_min, y_max) sets strict (hard) lower and 
%   upper bounds on the output. The solver will strictly enforce these 
%   limits. This is suitable for absolute physical boundaries, but may cause 
%   the solver to crash (go infeasible) if a disturbance pushes the system too far.
%
%   mpc = INIT_MPC_OUTPUT_CNSTR(mpc, y_min, y_max, y_min_slack_active, y_max_slack_active, qv_min, qv_max) 
%   allows you to define specific bounds as "soft" constraints. Soft constraints 
%   can be safely violated during massive disturbances to keep the solver running, 
%   while applying a customizable penalty to drive the state back within limits 
%   as quickly as possible.
%
%   INPUTS:
%       mpc                - CHRONOS MPC structure
%       y_min              - [ny x 1] Array of lower output limits (use [] if none).
%       y_max              - [ny x 1] Array of upper output limits (use [] if none).
%       qv_min             - (Optional) [ny x 1] or scalar. Penalty weight for violating 
%                            the y_min soft limits. Higher values mean stricter enforcement.
%       qv_max             - (Optional) [ny x 1] or scalar. Penalty weight for violating 
%                            the y_max soft limits. Higher values mean stricter enforcement.
%
%   OUTPUTS:
%       mpc                - Updated MPC structure. All necessary background math 
%                            (constraint gradients, Hessians, and slack variables) 
%                            are automatically assembled and added to the object.
%
%   USAGE TIPS:
%       - If qv_min or qv_max are not passed, the soft constraint penalty 
%         weight will default to the value stored in mpc.qv
%       - Passing a scalar to the slack or qv inputs will automatically apply 
%         that setting across all constrained outputs.
function mpc = init_mpc_output_cnstr(mpc,y_min,y_max,qv_min,qv_max)
arguments
    mpc
    y_min = []
    y_max = []
    qv_min = []
    qv_max = []
end

mpc.has_y_cnstr = 1;

% INPUT DIMENSION VALIDATION 
validate_column_vector(y_min, mpc.ny, 'y_min');
validate_column_vector(y_max, mpc.ny, 'y_max');
validate_column_vector(qv_min, mpc.ny, 'qv_min');
validate_column_vector(qv_max, mpc.ny, 'qv_max');

y_cnstr.use_k0 = mpc.y_use_k0;
y_cnstr.use_ter = mpc.y_use_ter;

% Expand scalars to full local vectors if needed
if isscalar(y_min), y_min = y_min * ones(mpc.ny, 1); end
if isscalar(y_max), y_max = y_max * ones(mpc.ny, 1); end

y_cnstr.min = y_min;
y_cnstr.max = y_max;

if mpc.y_use_k0
    if ~isempty(y_cnstr.min), y_cnstr.min_0 = y_min(mpc.y_rows_k0); end
    if ~isempty(y_cnstr.max), y_cnstr.max_0 = y_max(mpc.y_rows_k0); end
end

if mpc.y_use_ter
    if ~isempty(y_cnstr.min), y_cnstr.min_ter = y_min(mpc.y_rows_ter); end
    if ~isempty(y_cnstr.max), y_cnstr.max_ter = y_max(mpc.y_rows_ter); end
end


if ~isempty(y_cnstr.min)

    y_cnstr.min_limit = 1;
    
    if mpc.y_use_k0
        mpc.ng_k(1) = mpc.ng_k(1) + mpc.ny_0;
        mpc.nv_k(1) = mpc.nv_k(1) + mpc.ny_0;
    end
    mpc.ng_k(2) = mpc.ng_k(2) + mpc.ny;
    mpc.nv_k(2) = mpc.nv_k(2) + mpc.ny;
    if mpc.y_use_ter
        mpc.ng_k(3) = mpc.ng_k(3) + mpc.ny_ter;
        mpc.nv_k(3) = mpc.nv_k(3) + mpc.ny_ter;
    end

    y_cnstr.g_min_index_k = [];
    y_cnstr.v_min_index_k = [];

    % Initialize Penalty term for new slack variables
    if isempty(qv_min)
        % if qv isnt defined, it is not initialized until build_chronos_mpc(),
        % but we need to make space 
        if y_cnstr.use_k0, qv_min_0 = zeros(mpc.ny_0,1); end
        qv_min_k = zeros(mpc.ny,1);
        if y_cnstr.use_ter, qv_min_ter = zeros(mpc.ny_ter,1); end

    elseif isscalar(qv_min)
        if y_cnstr.use_k0, qv_min_0 = qv_min*ones(mpc.ny_0,1); end
        qv_min_k = qv_min*ones(mpc.ny,1);
        if y_cnstr.use_ter, qv_min_ter = qv_min*ones(mpc.ny_ter,1); end

    else % full vector is passed, pick elements for k=0 and k=N
        if y_cnstr.use_k0, qv_min_0 = qv_min(mpc.y_rows_k0); end
        qv_min_k = qv_min;
        if y_cnstr.use_ter, qv_min_ter = qv_min(mpc.y_rows_ter); end
    end

    if y_cnstr.use_k0, y_cnstr.qv_min_0 = qv_min_0; end
    y_cnstr.qv_min = qv_min_k;
    if y_cnstr.use_ter, y_cnstr.qv_min_ter = qv_min_ter; end

else
    y_cnstr.min_limit = 0;
end

if ~isempty(y_cnstr.max)

    y_cnstr.max_limit = 1;
    
    if mpc.y_use_k0
        mpc.ng_k(1) = mpc.ng_k(1) + mpc.ny_0;
        mpc.nv_k(1) = mpc.nv_k(1) + mpc.ny_0;
    end
    mpc.ng_k(2) = mpc.ng_k(2) + mpc.ny;
    mpc.nv_k(2) = mpc.nv_k(2) + mpc.ny;
    if mpc.y_use_ter
        mpc.ng_k(3) = mpc.ng_k(3) + mpc.ny_ter;
        mpc.nv_k(3) = mpc.nv_k(3) + mpc.ny_ter;
    end

    y_cnstr.g_max_index_k = [];
    y_cnstr.v_max_index_k = [];

    % Initialize Penalty term for new slack variables
    if isempty(qv_max)
        % if qv isnt defined, it is not initialized until build_chronos_mpc(),
        % but we need to make space 
        if y_cnstr.use_k0, qv_max_0 = zeros(mpc.ny_0,1); end
        qv_max_k = zeros(mpc.ny,1);
        if y_cnstr.use_ter, qv_max_ter = zeros(mpc.ny_ter,1); end

    elseif isscalar(qv_max)
        if y_cnstr.use_k0, qv_max_0 = qv_max*ones(mpc.ny_0,1); end
        qv_max_k = qv_max*ones(mpc.ny,1);
        if y_cnstr.use_ter, qv_max_ter = qv_max*ones(mpc.ny_ter,1); end

    else % full vector is passed, pick elements for k=0 and k=N
        if y_cnstr.use_k0, qv_max_0 = qv_max(mpc.y_rows_k0); end
        qv_max_k = qv_max;
        if y_cnstr.use_ter, qv_max_ter = qv_max(mpc.y_rows_ter); end
    end

    if y_cnstr.use_k0, y_cnstr.qv_max_0 = qv_max_0; end
    y_cnstr.qv_max = qv_max_k;
    if y_cnstr.use_ter, y_cnstr.qv_max_ter = qv_max_ter; end

else
    y_cnstr.max_limit = 0;
end

mpc.y_cnstr = y_cnstr;

end