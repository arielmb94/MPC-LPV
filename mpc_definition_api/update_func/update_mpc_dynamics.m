%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   mpc = update_mpc_sys_dynamics(mpc,A,B,Bd)
%
% Updates the matrices A, B and Bd of the internal MPC Discrete-Time system
% model:
% 
% x+ = A * x + B * u + Bd * d
% 
% and then recomputes the MPC equality constraints accordingly.
%
% Example uses:
%
%   - update only system matrix cost limits: 
%           mpc = update_mpc_sys_dynamics(mpc,A,[],[])
%   - update both system matrix and input matrox: 
%           mpc = update_mpc_sys_dynamics(mpc,[],B,[])
%   - update only input disturbance: 
%           mpc = update_mpc_sys_dynamics(mpc,[],[],Bd)
%
% In:
%   - mpc: CHRONOS mpc structure
%   - A (optional): nx x nx matrix, system matrix
%   - B (optional): nx x nu matrix, input matrix
%   - Bd (optional): nx x nd matrix, disturbance input matrix.
%
%   All arguments items which do not require to be updated can be passed as
%   an empty vector [].
%
% Out:
%   - mpc: updated CHRONOS mpc structure
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function mpc = update_mpc_dynamics(mpc,A,B,Bd)

if ~isempty(A)
    mpc.A = A;
    mpc.A_kkt(mpc.s_col,mpc.s_col) = mpc.A;
end

if ~isempty(B)
    mpc.B = B;
    mpc.B_kkt(mpc.s_col,:) = mpc.B;
end

if ~isempty(Bd)   
    mpc.Bd = Bd;
end

end