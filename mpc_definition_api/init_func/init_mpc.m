%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   mpc = init_mpc(N,N_ctr_hor)
%
% Initializes CHRONOS mpc structure fields and solver hyperparmeters
%
% In:
%   - N: MPC prediction horizon
%   - N_ctr_hor (optional):  prediction horizon for control actions. If not
%   specified it is set equal to N
%
% Out:
%   - mpc: initialized CRHONOS mpc structure
%
% Example Use:
%
%   - Same control and prediction horizons:
%               mpc = init_mpc(N)
%   - Different control and prediction horizons:  
%               mpc = init_mpc(N,N_ctr_hor)
%
% Hyperparameters (can be modified manually after initialization of the mpc
% structure):
%
%   - mpc.t: interior-point method tradeoff parameter between cost function
%   minimization vs constraint satisfaction. Large t values give preference
%   to minimization of the cost function. Small values for t will make the
%   solver prefer feasibility and constraint safety.
%
%   - mpc.Beta: reduction step for each iteration of the feasibility line
%   search. Beta must be less than 1 and greater than 0. Values close to 1
%   ensure a smoother optimization solution between multiple mpc
%   iterations at the cost of increased line search iterations.
%
%   - mpc.min_l: if the line search step fall below min_l the following
%   iteration of the interior-point method will be cancelled. Allows to
%   quit the interior-point method quicker when the optimal solution is
%   close to the constraints limits.
%
%   - mpc.eps: interior-point method precision.
%
%   - mpc.max_iter: maximum allowed iterations of the interior-point method
%   solver
%
%   - mpc.t_feas: exactly as mpc.t, applied for the step 0 feasibility
%   solver. The step 0 solver allows to find a feasibile starting point for
%   the interior-point method when providded the initial guess lies
%   outside of the feasible region.
%
%   - mpc.qfeas: cost term to penalize large deviation on the solution of
%   the step 0 solver from the provided initial guess
%
%   - mpc.v0_feas: initial value for step 0 solver slack variable
%
%   - mpc.feas_lambda: multiplier in case step 0 starting slack variable
%   value is set too low. Must be larger than 1.
%
%   - mpc.max_feas_iter: maximum number of step 0 solver iterations 
%   allowed. If max_feas_iter is violated it is assumed the problem is 
%   unfeasible. 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function mpc = init_mpc(N)
arguments
    N = 2
end

mpc.N = N;

mpc.Qe = [];
mpc.Rdu = [];
mpc.Ru = [];
mpc.ru = [];
mpc.A = [];
mpc.B = [];
mpc.Bd = [];
mpc.C = [];
mpc.D = [];
mpc.Dd = [];
mpc.Qz = [];
mpc.qz = [];
mpc.Cz = [];
mpc.Dz = [];
mpc.Ddz = [];
mpc.Dsuz = [];
mpc.Ch = [];
mpc.Dh = [];
mpc.Dsuh = [];
mpc.Ddh = [];
mpc.nx = 0;
mpc.nu = 0;
mpc.nd = 0;
mpc.ny = 0;
mpc.ndz = 0;
mpc.nz_0 = 0;
mpc.nz = 0;
mpc.nz_ter = 0;
mpc.ndh = 0;
mpc.nh_0 = 0;
mpc.nh = 0;
mpc.nh_ter = 0;
mpc.Nx = 0;
mpc.Nu = 0;
mpc.Nd = 0;

mpc.dyn_use_d = 0;
mpc.y_use_s = 0;
mpc.y_use_u = 0;
mpc.y_use_d = 0;
mpc.y_use_k0 = 0;
mpc.y_use_ter = 0;
mpc.z_use_s = 0;
mpc.z_use_u = 0;
mpc.z_use_su = 0;
mpc.z_use_d = 0;
mpc.z_use_k0 = 0;
mpc.z_use_ter = 0;

mpc.s = [];
mpc.s_ter = [];
mpc.u = [];
mpc.su = [];
mpc.du = [];
mpc.r = [];
mpc.y = [];
mpc.err = [];
mpc.d = [];
mpc.z = [];
mpc.dz = [];
mpc.h = [];
mpc.dh = [];
mpc.g = [];
mpc.v = [];
mpc.xN_ref = [];

mpc.tracking_cost = 0;
mpc.quad_control_cost = 0;
mpc.lin_control_cost = 0;
mpc.controlrate_cost = 0;
mpc.quad_custom_cost = 0;
mpc.lin_custom_cost = 0;
mpc.tracking_cost_index_k = [];
mpc.custom_cost_index_k = [];

mpc.update_tracking = 0;
mpc.update_customcost_quad = 0;
mpc.update_customcost_lin = 0;
mpc.recompute_cost_hess = 0;

mpc.s = [];
mpc.s_ter = [];
mpc.u = [];
mpc.du = [];
mpc.y = [];
mpc.h = [];
mpc.z = [];
mpc.nvar = 0;
mpc.t = 50;
mpc.Beta = 0.75;
mpc.min_l = 1e-6;
mpc.eps = 1e-4;
mpc.max_iter = 10;
mpc.ter_ingredients = 0;
mpc.xN_ref_is_y = 0;
mpc.P = [];
mpc.P2 = [];
mpc.K = [];
mpc.warm_starting = 0;

mpc.has_s_cnstr = 0;
mpc.has_u_cnstr = 0;
mpc.has_du = 0;
mpc.has_du_cnstr = 0;
mpc.has_y_cnstr = 0;
mpc.has_h_cnstr = 0;

mpc.s_cnstr = [];
mpc.u_cnstr = [];
mpc.du_cnstr = [];
mpc.y_cnstr = [];
mpc.h_cnstr = [];

mpc.ng_k = [0 0 0]; % inequalites per horizon step
mpc.nv_k = [0 0 0]; % soft inequalites per horizon step

mpc.qv = 50; % Slack variable penalty
mpc.Qv_fctr = 10;
mpc.v = [];
mpc.slack_epsilon = 1e-3;
mpc.eps_thknv = 1e-6;

end