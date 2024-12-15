function [mpc] = init_mpc_ter_ingredients_dlqr(mpc,Qx,Ru,...
    terminal_constraint,x_ref_is_y)


mpc.ter_ingredients = 1;
mpc.ter_constraint = terminal_constraint;
mpc.x_ref_is_y = x_ref_is_y;

[K,P] = dlqr(mpc.A,mpc.B,Qx,Ru);

mpc.P = P;

mpc.hessTerminalCost = zeros(mpc.Nx+mpc.Nu,mpc.Nx+mpc.Nu);
mpc.hessTerminalCost(end-mpc.nx+1: end,end-mpc.nx+1 : end) = P;

end