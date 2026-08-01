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

    len_A = size(A,3);
    if len_A < mpc.N
        mpc.A(:,:,:) = fill_mat(mpc.A, A, 1);
    else
        mpc.A(:,:,:) = A(:,:,1:mpc.N);
    end

    mpc.A_kkt(mpc.s_col,mpc.s_col,:) = mpc.A(:,:,2:mpc.N);
end

if ~isempty(B)

    len_Bd = size(B,3);
    if len_Bd < mpc.N
        mpc.B(:,:,:) = fill_mat(mpc.B, B, 1);
    else
        mpc.B(:,:,:) = B(:,:,1:mpc.N);
    end

    mpc.B_kkt_0(mpc.s_col,:) = mpc.B(:,:,1);
    mpc.B_kkt(mpc.s_col,:,:) = mpc.B(:,:,2:mpc.N);
end

if ~isempty(Bd)   

    len_Bd = size(Bd,3);
    if len_Bd < mpc.N
        mpc.Bd(:,:,:) = fill_mat(mpc.Bd, Bd, 1);
    else
        mpc.Bd(:,:,:) = Bd(:,:,1:mpc.N);
    end
end

end