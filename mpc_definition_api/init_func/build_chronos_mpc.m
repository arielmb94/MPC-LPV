% INIT_INITIAL_GUESS Generates a feasible warm-start primal vector for the MPC solver.
%
%   x0 = INIT_INITIAL_GUESS(mpc, s_prev, u_prev) calculates a strictly feasible 
%   initial guess for the primal optimization vector using a system rollout.
%   It mathematically guarantees that the initial guess respects input limits, 
%   rate limits, and hard state constraints, preventing solver crashes.
%
%   x0 = INIT_INITIAL_GUESS(mpc, s_prev, u_prev, x_ref, d_in, dh_in) allows the 
%   inclusion of reference trajectories and measured disturbances.
%
%   HOW IT WORKS:
%   The function simulates the system dynamics across the prediction horizon:
%     1. Control Strategy: If an LQR gain (mpc.K) and reference (x_ref) exist, 
%        it uses feedback (u = K*(x_ref - x)). Otherwise, it holds u_prev constant.
%     2. Input Clipping: Control actions are strictly bounded by absolute (mpc.u_cnstr)
%        and rate (mpc.du_cnstr) limits.
%     3. State Clamping: If hard state constraints (mpc.s_cnstr) are defined (without 
%        slacks), the state is clamped to remain strictly feasible, intentionally 
%        breaking the equality constraint to prevent an Interior Point method crash.
%     4. Soft Constraints: States with enabled slacks are allowed to violate bounds. 
%        The slack variables are then automatically sized to absorb the violation.
%
%   INPUTS:
%       mpc    - CHRONOS MPC structure.
%       s_prev - [nx x 1] Current measured state vector.
%       u_prev - [nu x 1] Last applied control input.
%       x_ref  - [nx x 1] (Optional) State reference. Defaults to [].
%       d_in   - [nd x 1] (Optional) Measured disturbances.
%       dh_in  - [ndh x 1] (Optional) Measured disturbance vector for custom constraints.
%
%   OUTPUTS:
%       x0     - Primal vector [u_0; x_1; u_1; ... x_Nc; ... x_N; v] 
function mpc = build_chronos_mpc(mpc,s_prev,u_prev,d_in,x_ref)
arguments
    mpc
    s_prev
    u_prev
    d_in = []
    x_ref = []
end

    % get index maps for optimization variables
    mpc = build_index(mpc);
    
    % init equality constraints
    mpc = genEqualities(mpc);

    % init soft slacks cost
    mpc = set_soft_cost_qv(mpc);

    % init costs
    mpc = init_costs(mpc);

    % compute primal variables vector
    len_d_in = size(d_in,2);
    if ~isempty(d_in) && len_d_in< mpc.N
        mpc.d(:,:) = fill_vec(mpc.d,d_in,1);
    else
        mpc.d(:,:) = d_in;
    end 

    x0 = rollstates(mpc,s_prev,u_prev,x_ref,mpc.d);

    x0(mpc.slack_index) = 1/mpc.t;
    mpc.x0 = x0;
    
end

function x0 = rollstates(mpc,s_prev,u_prev,x_ref,d_in)

x0 = zeros(mpc.n,1);

x_k = s_prev;
u_k_prev = u_prev;

for k = 1 : mpc.N
    % Compute Raw Control Input
    if ~isempty(mpc.K) && ~isempty(x_ref)
        % Option 2: Terminal ingredients exist and user pass x_ref
        u_raw = mpc.K * (x_ref - x_k);
    else
        % Option 1: Constant previous input
        u_raw = u_k_prev;
    end

    u_k = u_raw;

    % 2. Clip for Rate Constraints (Delta u)
    if mpc.has_du_cnstr
        % Check minimum rate limit
        if mpc.du_cnstr.min_limit
            du_min_strict = mpc.du_cnstr.min + mpc.slack_epsilon;
            u_k = max(u_k_prev + du_min_strict, u_k);
        end
        % Check maximum rate limit
        if mpc.du_cnstr.max_limit
            du_max_strict = mpc.du_cnstr.max - mpc.slack_epsilon;
            u_k = min(u_k_prev + du_max_strict, u_k);
        end
    end

    % 3. Clip for Absolute Constraints (u)
    if mpc.has_u_cnstr
        % Check minimum absolute limit
        if mpc.u_cnstr.min_limit
            u_min_strict = mpc.u_cnstr.min + mpc.slack_epsilon;
            u_k = max(u_min_strict, u_k);
        end
        % Check maximum absolute limit
        if mpc.u_cnstr.max_limit
            u_max_strict = mpc.u_cnstr.max - mpc.slack_epsilon;
            u_k = min(u_max_strict, u_k);
        end
    end

    % 4. Propagate Dynamics
    % x_{k+1} = A*x_k + B*u_k + D*d_k
    x_next = mpc.A * x_k + mpc.B * u_k;
    if ~isempty(mpc.Bd) && ~isempty(d_in)
        x_next = x_next + mpc.Bd * d_in(:,k);
    end
    % clamp x: for safety net in case we are dealing with unstable
    % system
    if mpc.has_s_cnstr
        % Check minimum state limits
        if mpc.s_cnstr.min_limit
            s_min_strict = mpc.s_cnstr.min + mpc.slack_epsilon;
            x_next = max(s_min_strict, x_next);
        end

        % Check maximum state limits
        if mpc.s_cnstr.max_limit
            s_max_strict = mpc.s_cnstr.max - mpc.slack_epsilon;
            x_next = min(s_max_strict, x_next);

        end
    end

    % 5. Map to Primal Optimization Vector
    % Indexing math for [u0; x1; u1; x2; ...]
    x0(mpc.u_index_k(:,k)) = u_k;
    if mpc.has_du && any(mpc.su_index_k(:,k))
        x0(mpc.su_index_k(:,k)) = u_k_prev;
    end
    x0(mpc.s_index_k(:,k+1)) = x_next;

    % 6. Prepare for next step
    x_k = x_next;
    u_k_prev = u_k;
end

end