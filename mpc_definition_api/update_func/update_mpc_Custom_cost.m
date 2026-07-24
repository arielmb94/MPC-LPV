%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   mpc = update_mpc_Lin_Custom_cost(mpc,Cz,Dz,Ddz,Qz,qz)
%
% Updates the parameters for the user-defined cost term and re-computes the
% MPC gradients and Hessians accordingly.
%
% The function can be used to update the model of the user-defined signal 
% z:
%
%   z = Cz * x + Dz * u + Ddz * dz
%
% by updating the matrices Cz, Dz and Ddz.
%
% The function can also be called to update the weight values for the
% quadratic Ru and linear ru penalty terms of the MPC cost functions:
%
%   J += z'*Qz*z + qz*z
%
% Example uses:
%
%   - update only quadratic cost penalty: 
%           mpc = update_mpc_Lin_Custom_cost(mpc,[],[],[],Qz)
%   - update only user-defined signal z model: 
%           mpc = update_mpc_Lin_Custom_cost(mpc,Cz,Dz,Ddz)
%   - update only input feedthrough Dz matrix and linear cost: 
%           mpc = update_mpc_Lin_Custom_cost(mpc,[],Di,[],[],qz)
%
% In:
%   - mpc: CHRONOS mpc structure
%   - Cz (optional): nz x nx matrix, state output matrix
%   - Dz (optional): nz x nu matrix, input feedtrhough output matrix
%   - Ddz (optional): nz x ndz matrix, disturbance feedtrhough output 
%   matrix
%   - Qz (optional): nz x nz square matrix, weights for the quadratic
%   penalty term on the user defined signal z
%   - qz (optional): nz column vector, weights for the linear penalty term
%   on the user defined signal z.
%
%   All arguments items which do not require to be updated can be passed as
%   an empty vector [].
%
% Out:
%   - mpc: updated CHRONOS mpc structure
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function mpc = update_mpc_Custom_cost(mpc,Qz,qz)

if ~isempty(Qz)
    mpc.update_customcost_quad = 1;
    mpc.Qz(:,:) = Qz;
    if mpc.z_use_k0, mpc.Qz_0(:,:) = Qz(mpc.z_rows_k0,mpc.z_rows_k0); end
    if mpc.z_use_ter, mpc.Qz_ter(:,:) = Qz(mpc.z_rows_ter,mpc.z_rows_ter); end
end
if ~isempty(qz)
    mpc.update_customcost_lin = 1;
    mpc.qz(:) = qz;
    if mpc.z_use_k0, mpc.qz_0(:) = qz(mpc.z_rows_k0); end
    if mpc.z_use_ter, mpc.qz_ter(:) = qz(mpc.z_rows_ter); end
end
    
end