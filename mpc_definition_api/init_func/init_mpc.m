% INIT_MPC Create the CHRONOS MPC structure.
%
%   mpc = INIT_MPC(N) creates an empty MPC structure with prediction
%   horizon N.
%
%   mpc = INIT_MPC() uses the default horizon N = 10.
%
%   Call this function first. Then define the dynamics, optionally replace
%   the default tracking output, add costs and constraints, and finish with
%   BUILD_CHRONOS_MPC.
%
%   Input:
%     N       - Prediction horizon. Default: 10.
%
%   Output:
%     mpc     - New CHRONOS MPC structure, ready for model, cost, and
%               constraint initialization.
%
%   Example:
%
%       mpc = init_mpc(20);
%       mpc = init_mpc_dynamics(mpc, A, B, []);
%       mpc = init_mpc_output(mpc, C, [], []);
function mpc = init_mpc(N)
arguments
    N = 10
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
mpc.v = [];
mpc.slack_epsilon = 1e-3;
mpc.eps_thknv = 1e-6;

end
