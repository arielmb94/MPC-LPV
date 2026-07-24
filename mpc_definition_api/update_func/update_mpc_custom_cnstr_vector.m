%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%    mpc = update_mpc_lin_custom_cnstr(mpc,Ch,Dh,Ddh,h_min,h_max)
%
% Updates the parameters for the user-defined custom constraint and
% re-computes the MPC gradients and Hessians accordingly.
% 
% The function can be used to update the model of the user-defined signal 
% h:
%
%   h = Ch * x + Dh * u + Ddh * dh
%
% by updating the matrices Ch, Dh and Ddh.
%
% The function can also be called to update the constraint limits:
%
%   h_min <= h <= h_max
%
% Example uses:
%
%   - update only constraint limits: 
%           mpc = init_mpc_delta_u_cnstr(mpc,[],[],[],h_min,h_max)
%   - update only user-defined signal h model: 
%           mpc = init_mpc_delta_u_cnstr(mpc,Ch,Dh,Ddh)
%   - update only input feedthrough Dh matrix : 
%           mpc = init_mpc_delta_u_cnstr(mpc,[],Dh,[])
%
% In:
%   - mpc: CHRONOS mpc structure
%   - Ch (optional): nh x nx matrix, state output matrix
%   - Dh (optional): nh x nu matrix, input feedtrhough matrix
%   - Ddh (optional): nh x ndh matrix, disturbance feedtrhough matrix
%   - h_min (optional): nh column vector, lower bound constraint values
%   on the user defined signal h
%   - h_max (optional): nh column vector, upper bound constraint values 
%   on the user defined signal h
%
%   All arguments items which do not require to be updated can be passed as
%   an empty vector [].
%
% Out:
%   - mpc: updated CHRONOS mpc structure
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function mpc = update_mpc_custom_cnstr_vector(mpc,Ch,Dh,Dsuh,Ddh)

if ~isempty(Ch)

    mpc.Ch(:,:) = Ch;

    if mpc.h_cnstr.use_k0, mpc.Ch_0(:,:) = Ch(mpc.h_cnstr.rows_k0,:); end
    if mpc.h_cnstr.use_ter, mpc.Ch_ter(:,:) = Ch(mpc.h_cnstr.rows_ter,:); end

    if mpc.h_cnstr.min_limit
        for k = 1:mpc.N-1
            mpc.Ai_k(mpc.h_cnstr.min_ineqRow_k,mpc.s_col,k) = -mpc.Ch;
        end
        if mpc.h_cnstr.use_ter
            mpc.Ai_ter(mpc.h_cnstr.min_ineqRow_ter,mpc.s_col) = -mpc.Ch_ter;
        end
    end
    if mpc.h_cnstr.max_limit
        for k = 1:mpc.N-1
            mpc.Ai_k(mpc.h_cnstr.max_ineqRow_k,mpc.s_col,k) = mpc.Ch;
        end
        if mpc.h_cnstr.use_ter
            mpc.Ai_ter(mpc.h_cnstr.max_ineqRow_ter,mpc.s_col) = mpc.Ch_ter;
        end
    end

end

if ~isempty(Dh)
    mpc.Dh(:,:) = Dh;
    % if there is D it means there is k0
    mpc.Dh_0(:,:) = Dh(mpc.h_cnstr.rows_k0,:); 

    if mpc.h_cnstr.min_limit
        mpc.Ai_0(mpc.h_cnstr.min_ineqRow_0,:) = -mpc.Dh_0;
        for k = 1:mpc.N-1
            mpc.Ai_k(mpc.h_cnstr.min_ineqRow_k,mpc.u_col,k) = -mpc.Dh;
        end 
    end
    if mpc.h_cnstr.max_limit
        mpc.Ai_0(mpc.h_cnstr.max_ineqRow_0,:) = mpc.Dh_0;
        for k = 1:mpc.N-1
            mpc.Ai_k(mpc.h_cnstr.max_ineqRow_k,mpc.u_col,k) = mpc.Dh;
        end 
    end

end

if ~isempty(Dsuh)
    mpc.Dsuh(:,:) = Dsuh;
    if mpc.h_cnstr.use_k0, mpc.Dsuh_0(:,:) = Dsuh(mpc.h_cnstr.rows_k0,:); end

    if mpc.h_cnstr.min_limit
        for k = 1:mpc.N-1
            mpc.Ai_k(mpc.h_cnstr.min_ineqRow_k,mpc.su_col,k) = -mpc.Dsuh;
        end 
    end
    if mpc.h_cnstr.max_limit
        for k = 1:mpc.N-1
            mpc.Ai_k(mpc.h_cnstr.max_ineqRow_k,mpc.su_col,k) = mpc.Dsuh;
        end 
    end
end

if ~isempty(Ddh)   
    mpc.Ddh(:,:) = Ddh;
    if mpc.h_cnstr.use_k0, mpc.Ddh_0(:,:) = Ddh(mpc.h_cnstr.rows_k0,:); end
end

end