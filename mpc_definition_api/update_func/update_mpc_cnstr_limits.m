%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   mpc = update_mpc_u_cnstr(mpc,u_min,u_max)
%
% Modifies the constraints limits on the control action
%
% In:
%   - mpc: CHRONOS mpc structure
%   - u_min (optional): nu column vector, lower bound constraint values on
%   the control action
%   - u_max (optional): nu column vector, upper bound constraint values on
%   the control action 
%
% Out:
%   - mpc: updated CHRONOS mpc structure
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cnstr = update_mpc_cnstr_limits(cnstr,min,max)

if ~isempty(min)    
    cnstr.min(:) = min;
    if cnstr.use_k0
        cnstr.min_0(:) = min(cnstr.rows_k0);
    end
    if cnstr.use_ter
        cnstr.min_ter(:) = min(cnstr.rows_ter);
    end
end

if ~isempty(max)    
    cnstr.max(:) = max;
    if cnstr.use_k0
        cnstr.max_0(:) = max(cnstr.rows_k0);
    end
    if cnstr.use_ter
        cnstr.max_ter(:) = max(cnstr.rows_ter);
    end
end

end