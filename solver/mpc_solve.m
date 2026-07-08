%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
% Solve the current iteration of the MPC problem.
%
% In:
%   - mpc: CHRONOS mpc structure.
%   - x0: Nx+Nu column vector, initial guess solution for CHRONOS interior
%   point iterative solver
%   - s_prev: nx column vector, last measured or estimated system state
%   value
%   - u_prev: nu column vector, control action applied to the system on the
%   previous sampling time
%   - r (optional): tracking reference for the MPC. It can be a ny column
%   vector (the same reference applies for the full prediction horizon) or
%   can be a Ny column vector (the user passes a unique reference for each
%   step of the prediction horizon). If not used, the user must pass an
%   empty vector [].
%   - d (optional): disturbance input to the system dynamics and to the
%   output signal y models. It can be a nd column vector (the same
%   disturbance applies for the full prediction horizon) or can be an Nd
%   column vector (the user passes a unique disturbance for each step of
%   the prediction horizon). If not used, the user must pass an empty
%   vector [].
%   - x_ref (optional): nx column vector, reference for the terminal state
%   xN of the prediction horizon. Required whenever the MPC problem
%   contains terminal ingredients. If not used, the user must pass an empty
%   vector [].
%   - dz (optional): disturbance input to the user defined signal model z
%   for custom cost functions. It can be a ndz column vector (the same
%   disturbance applies for the full prediction horizon) or can be an Ndz
%   column vector (the user passes a unique disturbance for each step of
%   the prediction horizon). If not used, the user must pass an empty
%   vector [].
%   - dh (optional): disturbance input to the user defined signal model h
%   for custom constraints. It can be a ndh column vector (the same
%   disturbance applies for the full prediction horizon) or can be an Ndh
%   column vector (the user passes a unique disturbance for each step of
%   the prediction horizon). If not used, the user must pass an empty
%   vector [].
%
% Out:
%   - u0: nu column vector, first step of the control action sequence
%   computed as solution to the MPC problem
%   - x0: Nx+Nu column vector, optimization variables solution vector to
%   the MPC problem
%   - iter: number of iterations required for the MPC optimization problem
%   - iter_feas: number of iterations required for the step 0 feasibility
%   starting point finder
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [u0,iter,mpc] = mpc_solve(mpc,s_prev,u_prev, ...
                                   r_in,xN_ref_in,...
                                   d_in,dz_in,dh_in)

% number of variables
n = mpc.n;
% number of equality constraints
n_eq = size(mpc.Aeq,1);

x0 = mpc.x0;
% handle input vector sizes
len_r_in = size(r_in,2);
if ~isempty(r_in) && len_r_in < mpc.N-1
    mpc.r(:,:) = fill_vec(mpc.r,r_in,1);
else
    mpc.r(:,:) = r_in;
end

len_d_in = size(d_in,2);
if ~isempty(d_in) && len_d_in < mpc.N
    mpc.d(:,:) = fill_vec(mpc.d,d_in,1);
else
    mpc.d(:,:) = d_in;
end

len_dz_in = size(dz_in,2);
if ~isempty(dz_in) && len_dz_in < mpc.Nz
    mpc.dz(:,:) = fill_vec(mpc.dz,dz_in,1);
else
    mpc.dz(:,:) = dz_in;
end

len_dh_in = size(dh_in,2);
if ~isempty(dh_in) && len_dh_in < mpc.Nh
    mpc.dh(:,:) = fill_vec(mpc.dh,dh_in,1);
else
    mpc.dh(:,:) = dh_in;
end

if mpc.ter_ingredients
    if mpc.xN_ref_is_y && isempty(xN_ref_in)
        mpc.xN_ref(:) = mpc.r(:,mpc.N-1);
    else
        mpc.xN_ref(:) = xN_ref_in;
    end
end

% update b matrix from equality condition
mpc = update_mpc_beq(mpc,s_prev,u_prev);

% Recompute hessian if cost terms have been updated
if mpc.recompute_cost_hess
    mpc = update_mpc_f0_hess(mpc);
end

% Set Newton solver condition at start
continue_Newton = true;
iter = 0;

mpc = get_mpc_variables(mpc,x0,s_prev,u_prev);

opts.SYM = true;
lambda2 = 1;

while mpc.eps <= lambda2*0.5 && continue_Newton && iter < mpc.max_iter

    % Compute gradient:

    % 1. Compute gradient/Hessian of box inequalities at x0:
    % init inequalities gradient vector
    grad_fi_Ind = zeros(n,1);
    grad_fi_Ind(mpc.slack_index) = -1./(mpc.slacks-mpc.slack_epsilon);
    % init inequalities hessian vector
    hess_fi_Ind = zeros(n,1);
    hess_fi_Ind(mpc.slack_index) = 1./(mpc.slacks-mpc.slack_epsilon).^2;

    mpc = grad_f0_MPC(mpc);

    mpc = equality_residuals(mpc);

    mpc = reduced_KKT_elements(mpc);

    [delta_u,delta_se,mu] = riccati_KKT(mpc,mpc.Q_k,mpc.R_k,mpc.Y_k,...  
                            mpc.ru_hat_k,mpc.rse_hat_k,mpc.rp)

    % 4. Compute gradient at x0 : grad(J) = t*grad(f0)+grad(Phi)
    grad_J_x0 = mpc.t*grad_f0+grad_fi_Ind;

    % 3. Compute Hessian of f(x0,t):
    hess_J_x0 = mpc.t*mpc.hessCost+mpc.eps_thknv*eye(n);
    for k = 1:length(mpc.slack_index)
        i = mpc.slack_index(k);
        hess_J_x0(i,i) = hess_J_x0(i,i) + hess_fi_Ind(i);
    end

    % solve KKT system
    %KKT = [hess_J_x0 mpc.Aeq';mpc.Aeq zeros(n_eq)];

   [delta_var,delta_g,delta_v] = reduced_KKT(mpc,x0,grad_f0,opts);
   delta_x_prim = zeros(n,1);
   delta_x_prim(mpc.variables_index) = delta_var;
   delta_x_prim(mpc.g_index) = delta_g;
   delta_x_prim(mpc.v_index) = delta_v;

%     delta_x = - linsolve(KKT,[grad_J_x0;mpc.Aeq*x0-mpc.beq],opts);
%     delta_x_prim = delta_x(1:n);

    % compute lambda^2
    lambda2 = -grad_J_x0'*delta_x_prim;

    % Feasibility line search
    l = 1;
    xhat = x0+l*delta_x_prim;

    feas = all(xhat(mpc.slack_index)>mpc.slack_epsilon);

    if feas
        x0 = xhat;
    else
        while ~feas
            l = l*mpc.Beta;

            xhat = x0+l*delta_x_prim;

            feas = all(xhat(mpc.slack_index)>mpc.slack_epsilon);
        end
        x0 = xhat;
        if l<mpc.min_l
            continue_Newton = false;
        end
    end
    mpc = get_mpc_variables(mpc,x0,s_prev,u_prev);
    iter = iter+1;
end

u0 = mpc.u(:,1);
mpc.x0(:) = x0;


end