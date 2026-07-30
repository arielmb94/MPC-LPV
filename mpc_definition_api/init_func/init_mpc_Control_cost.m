%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   mpc = init_mpc_Control_cost(mpc,Ru,ru)
%
% Adds quadratic and linear penalties on the control action:
%
%   J += u'*Ru*u + ru*u
%
% The function can be called to initialize either the quadratic term,
% either the linear term or both.
%
% Example uses:
%
%   - only linear term: mpc = init_mpc_Control_cost(mpc,[],ru)
%   - only quadratic term: mpc = init_mpc_Control_cost(mpc,Ru)
%   - both quadratic and linear terms: mpc = init_mpc_Control_cost(mpc,Ru,ru)
%
% In:
%   - mpc: CHRONOS mpc structure.
%   - Ru (optional): nu x nu square matrix, weights for the quadratic
%   penalty term on the control action.
%   - ru (optional): nu column vector, weights for the linear penalty term
%   on the control action. IMPORTANT: Use linear penalties only in the case
%   that the control action takes positive values only.
%
% Out:
%   - mpc: updated CHRONOS mpc structure
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function mpc = init_mpc_Control_cost(mpc,Ru,ru)
arguments
    mpc
    Ru = []
    ru = [];
end

if any(Ru(:))
    mpc.quad_control_cost = 1;

    mpc.Ru = zeros(mpc.nu,mpc.nu,mpc.N);
    len_Ru = size(Ru,3);
    if len_Ru < mpc.N
        mpc.Ru = fill_mat(mpc.Ru, Ru, 1);
    else
        mpc.Ru = Ru(:,:,1:mpc.N);
    end
end

if any(ru(:))
    mpc.lin_control_cost = 1;

    mpc.ru = zeros(mpc.nu,mpc.N);
    len_ru = size(ru,2);
    if len_ru < mpc.N
        mpc.ru = fill_vec(mpc.ru, ru, 1);
    else
        mpc.ru = ru(:,1:mpc.N);
    end    
end

end