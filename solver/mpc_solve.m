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
function [u0,iter,mpc] = mpc_solve(mpc,s_prev,u_prev,r_in,xN_ref_in,...
                                   d_in,dz_in,dh_in)

x0 = mpc.x0;

% handle input vector sizes
len_r_in = size(r_in,2);
if len_r_in 
    if len_r_in < mpc.N
        mpc.r(:,:) = fill_vec(mpc.r,r_in,1);
        if mpc.y_use_k0, mpc.r_0(:) = mpc.r(:,1); end
        if mpc.y_use_ter, mpc.r_ter(:) = mpc.r(:,mpc.N-1); end
    else 
        mpc.r(:,:) = r_in(:,1:mpc.N-1);
        if mpc.y_use_k0, mpc.r_0(:) = mpc.r(:,1); end
        if mpc.y_use_ter, mpc.r_ter(:) = r_in(:,mpc.N); end
    end
end

len_d_in = size(d_in,2);
if ~isempty(d_in) && len_d_in < mpc.N
    mpc.d(:,:) = fill_vec(mpc.d,d_in,1);
else
    mpc.d(:,:) = d_in;
end

len_dz_in = size(dz_in,2);
if ~isempty(dz_in) && len_dz_in < mpc.N
    mpc.dz(:,:) = fill_vec(mpc.dz,dz_in,1);
else
    mpc.dz(:,:) = dz_in;
end

len_dh_in = size(dh_in,2);
if ~isempty(dh_in) && len_dh_in < mpc.N
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

% Recompute gradient/hessian if cost terms have been updated
if mpc.update_tracking
    mpc = update_tracking_cost(mpc);
end
if mpc.update_customcost_quad
    mpc = update_custom_cost_quad(mpc);
end
if mpc.update_customcost_lin
    mpc = update_custom_cost_lin(mpc);
end
if mpc.recompute_cost_hess
    mpc = update_mpc_f0_hess(mpc);
end

% update dynamics equality RHS
mpc = update_mpc_beq(mpc,s_prev);

% Set Newton solver condition at start
continue_Newton = true;
iter = 0;

mpc = get_mpc_variables(mpc,x0,s_prev,u_prev);

lambda2 = 1;

while mpc.eps <= lambda2*0.5 && continue_Newton && iter < mpc.max_iter

    [mpc.grad_u_f0_0,...
     mpc.grad_se_f0_k,mpc.grad_u_f0_k,...
     mpc.grad_se_f0_ter] = grad_f0_MPC(mpc,mpc.grad_u_f0_0,mpc.grad_se_f0_k,...
                                       mpc.grad_u_f0_k,mpc.grad_se_f0_ter,...
                                       mpc.s_col,mpc.su_col,mpc.N,mpc.t);

    [mpc.rp_0,mpc.rp_k] = equality_residuals(mpc.rp_0,mpc.rp_k,mpc.A,mpc.B,...
                            mpc.beq_0,mpc.beq_k,mpc.u,mpc.s,...
                            mpc.has_du,mpc.su,mpc.s_col,mpc.su_col,mpc.N);

    if any(mpc.ng_k)
        [mpc.ri_0,mpc.ri_k,mpc.ri_ter] = inequality_residuals(mpc,mpc.N,...
                                        mpc.ri_0,mpc.ri_k,mpc.ri_ter,...
                                        mpc.g_0,mpc.g_k,mpc.g_ter,...
                                        mpc.v_0,mpc.v_k,mpc.v_ter,...
                                        mpc.ng_k,mpc.nv_k,mpc.v_rows_0,mpc.v_rows_k,...
                                        mpc.s_cnstr,mpc.u_cnstr,mpc.du_cnstr,...
                                        mpc.y_cnstr,mpc.h_cnstr);
    end

    [mpc.ru_0,mpc.rse_k,mpc.ru_k,mpc.rse_ter,...
     mpc.R_0,mpc.Q_k,mpc.R_k,mpc.Y_k,mpc.Q_ter,...
     mpc.g2_0,mpc.g2_k,mpc.g2_ter,mpc.v2_0,mpc.v2_k,mpc.v2_ter,...
     mpc.rv_v2_0,mpc.rv_v2_k,mpc.rv_v2_ter,...
     mpc.ri_0,mpc.ri_k,mpc.ri_ter,...
     mpc.iS_0,mpc.iS_k,mpc.iS_ter,...
     mpc.iS_ri_hat_0,mpc.iS_ri_hat_k,mpc.iS_ri_hat_ter] = reduced_KKT_elements(mpc,...
                        mpc.t,mpc.N,mpc.nu,mpc.nx,mpc.R_0,mpc.Q_k,mpc.R_k,mpc.Y_k,mpc.Q_ter,...
                        mpc.R_f0_0,mpc.Q_f0_k,mpc.R_f0_k,mpc.Y_f0_k,mpc.Q_f0_ter,...
                        mpc.ru_0,mpc.rse_k,mpc.ru_k,mpc.rse_ter,...
                        mpc.grad_u_f0_0,mpc.grad_se_f0_k,mpc.grad_u_f0_k,mpc.grad_se_f0_ter,...
                        mpc.u_cnstr,mpc.du_cnstr,mpc.s_cnstr,mpc.s_col,mpc.su_col,...
                        mpc.ng_k,mpc.nv_k,mpc.g_0,mpc.g_k,mpc.g_ter,...
                        mpc.v_0,mpc.v_k,mpc.v_ter,mpc.g2_0,mpc.g2_k,mpc.g2_ter,...
                        mpc.v2_0,mpc.v2_k,mpc.v2_ter,mpc.rv_v2_0,mpc.rv_v2_k,mpc.rv_v2_ter,...
                        mpc.ri_0,mpc.ri_k,mpc.ri_ter,mpc.grad_qv_0,mpc.grad_qv_k,...
                        mpc.grad_qv_ter,mpc.v_rows_0,mpc.v_rows_k,...
                        mpc.iS_0,mpc.iS_k,mpc.iS_ter,...
                        mpc.iS_ri_hat_0,mpc.iS_ri_hat_k,mpc.iS_ri_hat_ter,...
                        mpc.y_cnstr,mpc.D_0,mpc.C,mpc.D,mpc.C_ter,...
                        mpc.h_cnstr,mpc.Dh_0,mpc.Ch,mpc.Dsuh,mpc.Dh,mpc.Ch_ter);

    [mpc.delta_u,mpc.delta_se,...
     mpc.Q_hat,mpc.R_hat_0,mpc.R_hat,mpc.Y_hat,...
     mpc.rs_hat,mpc.ru_hat,mpc.ru_hat_0,...
     mpc.rp_hat_0,mpc.rp_hat] = riccati_KKT(mpc.N,...
                        mpc.R_0,mpc.Q_k,mpc.R_k,mpc.Y_k,mpc.Q_ter,...
                        mpc.Q_hat,mpc.R_hat_0,mpc.R_hat,mpc.Y_hat,...
                        mpc.ru_0,mpc.rse_k,mpc.ru_k,mpc.rse_ter,...
                        mpc.rs_hat,mpc.ru_hat,mpc.ru_hat_0,...
                        mpc.rp_0,mpc.rp_k,mpc.rp_hat_0,mpc.rp_hat,...
                        mpc.A_kkt,mpc.B_kkt,mpc.B_kkt_0,...
                        mpc.delta_u,mpc.delta_se,...
                        mpc.eps_thknv,mpc.nu,mpc.nse);

    if any(mpc.ng_k)
        [mpc.delta_g_0,mpc.delta_g_k,mpc.delta_g_ter,...
         mpc.delta_v_0,mpc.delta_v_k,mpc.delta_v_ter] = recover_slacks(mpc,mpc.N,...
                    mpc.delta_g_0,mpc.delta_g_k,mpc.delta_g_ter,...
                    mpc.delta_v_0,mpc.delta_v_k,mpc.delta_v_ter,...
                    mpc.delta_se,mpc.delta_u,mpc.iS_ri_hat_0,mpc.iS_ri_hat_k,...
                    mpc.iS_ri_hat_ter,mpc.iS_0,mpc.iS_k,mpc.iS_ter,...
                    mpc.g_0,mpc.g_k,mpc.g_ter,mpc.g2_0,mpc.g2_k,mpc.g2_ter,...
                    mpc.v2_0,mpc.v2_k,mpc.v2_ter,...
                    mpc.rv_v2_0,mpc.rv_v2_k,mpc.rv_v2_ter,...
                    mpc.s_cnstr,mpc.u_cnstr,mpc.du_cnstr,...
                    mpc.nu,mpc.nx,mpc.su_col,...
                    mpc.y_cnstr,mpc.D_0,mpc.C,mpc.D,mpc.C_ter,mpc.s_col,...
                    mpc.h_cnstr,mpc.Dh_0,mpc.Ch,mpc.Dsuh,mpc.Dh,mpc.Ch_ter);
    end

    [delta_x_prim,grad_J_x0] = stage2vec(mpc.delta_u,mpc.delta_se,mpc.delta_g_0,...
                                mpc.delta_g_k,mpc.delta_g_ter,...
                                mpc.delta_v_0,mpc.delta_v_k,mpc.delta_v_ter,...
                                mpc.g_0,mpc.g_k,mpc.g_ter,...
                                mpc.v_0,mpc.v_k,mpc.v_ter,...
                                mpc.t,mpc.grad_qv_0,mpc.grad_qv_k,mpc.grad_qv_ter,...
                                mpc.grad_u_f0_0,mpc.grad_u_f0_k,...
                                mpc.grad_se_f0_k,mpc.grad_se_f0_ter,...
                                mpc.u_index_k,mpc.se_index_k,...
                                mpc.g_index_0,mpc.g_index_k,mpc.g_index_ter,...
                                mpc.v_index_0,mpc.v_index_k,mpc.v_index_ter,...
                                mpc.ng_k,mpc.nv_k,mpc.n,mpc.N);

    % compute lambda^2
    lambda2 = -grad_J_x0*delta_x_prim;

    % Feasibility line search
    l = 1;
    xhat = x0+l*delta_x_prim;

    feas = all(xhat(mpc.g_index)>0) &&...
           all(xhat(mpc.v_index)>0);

    if feas
        x0 = xhat;
    else
        while ~feas
            l = l*mpc.Beta;

            xhat = x0+l*delta_x_prim;

            feas = all(xhat(mpc.g_index)>0) &&...
                   all(xhat(mpc.v_index)>0);
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
