%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   mpc = update_mpc_sys_output(mpc,C,D,Dd,Qe,y_min,y_max)
%
% Allows to update all parameters related to the tracking output signal y 
% in a single function:
%
%   - Upadate output signal model: y = C * x + D * u + Dd * d
%   - Upadate weight Qe on the tracking error penalty term: 
%   (r - y)' * Qe * (r - y)
%   - Upadate contraint Limits on feedback signal: y_min <= y <= y_max
%
% Example uses:
%
%   - update only the output feedback signal model: 
%           mpc = update_mpc_sys_output(mpc,C,D,Dd)
%   - update only the input feedtrhough matrix of the output signal model: 
%           mpc = update_mpc_sys_output(mpc,[],D,[])
%   - update the feedback output signal model and constraint limits: 
%           mpc = update_mpc_sys_output(mpc,C,D,Dd,[],y_min,y_max)
%   - update only the weight on the tracking error penalty term: 
%           mpc = update_mpc_sys_output(mpc,[],[],[],Qe)
%
% In:
%   - mpc: CHRONOS mpc structure
%   - C (optional): ny x nx matrix, system output matrix
%   - D (optional): ny x nu matrix, input feedtrhough matrix.
%   - Dd (optional): ny x nd matrix, disturbance feedtrhough matrix.
%   - Qe (optional): ny x ny square matrix, weights for the quadratic
%   penalty on the tracking error
%   - y_min (optional): ny column vector, lower bound constraint values on 
%   the tracking signal
%   - y_max (optional): ny column vector, upper bound constraint values on 
%   the tracking signal
%
%   All arguments items which do not require to be updated can be passed as
%   an empty vector [].
%
% Out:
%   - mpc: updated CHRONOS mpc structure
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function mpc = update_mpc_Custom_cost_vector(mpc,Cz,Dz,Dsuz,Ddz)

mpc.update_customcost_quad = mpc.quad_custom_cost;
mpc.update_customcost_lin = mpc.lin_custom_cost;
update_grad = 0;

if ~isempty(Cz)
    update_grad = 1;

    len_Cz = size(Cz,3);
    if len_Cz < mpc.N
        mpc.Cz(:,:,:) = fill_mat(mpc.Cz, Cz, 1);
        if mpc.z_use_ter, mpc.Cz_ter(:,:) = mpc.Cz(mpc.z_rows_ter,:,mpc.N-1); end
    else
        mpc.Cz(:,:,:) = Cz(:,:,1:mpc.N-1);
        if mpc.z_use_ter, mpc.Cz_ter(:,:) = Cz(mpc.z_rows_ter,:,mpc.N); end
    end
    if mpc.z_use_k0, mpc.Cz_0(:,:) = mpc.Cz(mpc.z_rows_k0,:,1); end

    if mpc.z_use_ter
        mpc.grad_z_ter(:,:) = mpc.Cz_ter';
    end
end

if ~isempty(Dz)
    update_grad = 1;

    len_Dz = size(Dz,3);
    if len_Dz < mpc.N-1
        mpc.Dz(:,:,:) = fill_mat(mpc.Dz, Dz, 1);
    else
        mpc.Dz(:,:,:) = Dz(:,:,1:mpc.N-1);
    end
    % if there is Dz it means there is k0
    mpc.Dz_0(:,:) = mpc.Dz(mpc.z_rows_k0,:,1);

    mpc.grad_z_0(:,:) = mpc.Dz_0';
end

if ~isempty(Dsuz)
    update_grad = 1;
    
    len_Dsuz = size(Dsuz,3);
    if len_Dsuz < mpc.N-1
        mpc.Dsuz(:,:,:) = fill_mat(mpc.Dsuz, Dsuz, 1);
    else
        mpc.Dsuz(:,:,:) = Dsuz(:,:,1:mpc.N-1);
    end
    if mpc.z_use_k0, mpc.Dsuz_0(:,:) = mpc.Dsuz(mpc.z_rows_k0,:,1); end
end

if ~isempty(Ddz)   

    len_Ddz = size(Ddz,3);
    if len_Ddz < mpc.N-1
        mpc.Ddz(:,:,:) = fill_mat(mpc.Ddz, Ddz, 1);
    else
        mpc.Ddz(:,:,:) = Ddz(:,:,1:mpc.N-1);
    end
    if mpc.z_use_k0, mpc.Ddz_0(:,:) = mpc.Ddz(mpc.z_rows_k0,:,1); end
end

if update_grad

    for k = 1:mpc.N-1
        if mpc.z_use_s && mpc.z_use_su && mpc.z_use_u
            mpc.grad_z(:,:,k) = [mpc.Cz(:,:,k)'; mpc.Dsuz(:,:,k)'; mpc.Dz(:,:,k)'];
        elseif mpc.z_use_s && mpc.z_use_su
            mpc.grad_z(:,:,k) = [mpc.Cz(:,:,k)'; mpc.Dsuz(:,:,k)'];
        elseif mpc.z_use_s && mpc.z_use_u
            mpc.grad_z(:,:,k) = [mpc.Cz(:,:,k)'; mpc.Dz(:,:,k)'];
        elseif mpc.z_use_su && mpc.z_use_u
            mpc.grad_z(:,:,k) = [mpc.Dsuz(:,:,k)'; mpc.Dz(:,:,k)'];
        elseif mpc.z_use_s
            mpc.grad_z(:,:,k) = mpc.Cz(:,:,k)';
        elseif mpc.z_use_su
            mpc.grad_z(:,:,k) = mpc.Dsuz(:,:,k)';
        elseif mpc.z_use_u
            mpc.grad_z(:,:,k) = mpc.Dz(:,:,k)';
        end
    end
end

end