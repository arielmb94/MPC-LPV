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
function mpc = update_mpc_output_vector(mpc,C,D,Dd)

mpc.update_tracking = 1;
update_grad = 0;

if ~isempty(C)
    update_grad = 1;
    mpc.C(:,:) = C;
    if mpc.y_use_k0, mpc.C_0(:,:) = C(mpc.y_rows_k0,:); end
    if mpc.y_use_ter
        mpc.C_ter(:,:) = C(mpc.y_rows_ter,:); 
        mpc.grad_err_ter(:,:) = -mpc.C_ter';
    end

    if mpc.has_y_cnstr
        if mpc.y_cnstr.min_limit
            for k = 1:mpc.N-1
                mpc.Ai_k(mpc.y_cnstr.min_ineqRow_k,mpc.s_col,k) = -mpc.C;
            end
            if mpc.y_cnstr.use_ter
                mpc.Ai_ter(mpc.y_cnstr.min_ineqRow_ter,mpc.s_col) = -mpc.C_ter;
            end
        end
        if mpc.y_cnstr.max_limit
            for k = 1:mpc.N-1
                mpc.Ai_k(mpc.y_cnstr.max_ineqRow_k,mpc.s_col,k) = mpc.C;
            end
            if mpc.y_cnstr.use_ter
                mpc.Ai_ter(mpc.y_cnstr.max_ineqRow_ter,mpc.s_col) = mpc.C_ter;
            end
        end
    end

end

if ~isempty(D)
    update_grad = 1;
    mpc.D(:,:) = D;
    % if there is D it means there is k0
    mpc.D_0(:,:) = D(mpc.y_rows_k0,:); 
    mpc.grad_err_0(:,:) = -mpc.D_0';

    if mpc.has_y_cnstr
        if mpc.y_cnstr.min_limit
            mpc.Ai_0(mpc.y_cnstr.min_ineqRow_0,:) = -mpc.D_0;
            for k = 1:mpc.N-1
                mpc.Ai_k(mpc.y_cnstr.min_ineqRow_k,mpc.u_col,k) = -mpc.D;
            end
        end
        if mpc.y_cnstr.max_limit
            mpc.Ai_0(mpc.y_cnstr.max_ineqRow_0,:) = mpc.D_0;
            for k = 1:mpc.N-1
                mpc.Ai_k(mpc.y_cnstr.max_ineqRow_k,mpc.u_col,k) = mpc.D;
            end
        end
    end

end

if ~isempty(Dd)   
    mpc.Dd(:,:) = Dd;
    if mpc.y_use_k0, mpc.Dd_0(:,:) = Dd(mpc.y_rows_k0,:); end
end

if update_grad
    if mpc.y_use_s && mpc.y_use_u
        mpc.grad_err(:,:) = [-mpc.C'; -mpc.D'];
    elseif mpc.y_use_s
        mpc.grad_err(:,:) = -mpc.C';
    elseif mpc.y_use_u
        mpc.grad_err(:,:) = -mpc.D';
    end
end

end