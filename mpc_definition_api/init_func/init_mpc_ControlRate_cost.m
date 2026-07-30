%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%    mpc = init_mpc_DiffControl_cost(mpc,Rdu)
%
% Adds penalty on the control variation between sampling instances:
%
%   J += delta_u' * Rdu * delta_u
%
% In:
%   - mpc: CHRONOS mpc structure
%   - Rdu: nu x nu square matrix, weights for the quadratic penalty term on
%   the control control variation between sampling instances
%
% Out:
%   - mpc: updated CHRONOS mpc structure
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function mpc = init_mpc_ControlRate_cost(mpc,Rdu)

mpc.has_du = 1;
mpc.controlrate_cost = 1;

mpc.Rdu = zeros(mpc.nu,mpc.nu,mpc.N);
len_Rdu = size(Rdu,3);
if len_Rdu < mpc.N
    mpc.Rdu = fill_mat(mpc.Rdu, Rdu, 1);
else
    mpc.Rdu = Rdu(:,:,1:mpc.N);
end

end