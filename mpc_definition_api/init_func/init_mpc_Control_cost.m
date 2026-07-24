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

mpc.Ru = Ru;
mpc.ru = ru;

if ~isempty(Ru)
    mpc.quad_control_cost = 1;
end
if ~isempty(ru)
    mpc.lin_control_cost = 1;
end

end