function [Aeq, beq] = appendGeneralizedConstraint(Aeq, beq, row, s_col, u_col, su_col, g_col, v_col, C_mat, D_mat, Ddu_mat, g_mat, v_mat, b_val)
% APPENDGENERALIZEDCONSTRAINT Inserts a unified linear constraint block.
% Form: C*s + D*u + Ddu*su + g_mat*g + v_mat*v = b_val

    % 1. State term (C * s_k)
    if ~isempty(s_col) && ~isempty(C_mat)
        Aeq(row, s_col) = C_mat;
    end
    
    % 2. Control Input term (D * u_k)
    if ~isempty(u_col) && ~isempty(D_mat)
        Aeq(row, u_col) = D_mat;
    end
    
    % 3. Rate/Past Input term (Ddu * su_k)
    if ~isempty(su_col) && ~isempty(Ddu_mat)
        Aeq(row, su_col) = Ddu_mat;
    end
    
    % 4. Slack translation term (g_mat * g_k)
    if ~isempty(g_col) && ~isempty(g_mat)
        Aeq(row, g_col) = g_mat;
    end
    
    % 5. Soft constraint relaxation term (v_mat * v_k)
    if ~isempty(v_col) && ~isempty(v_mat)
        Aeq(row, v_col) = v_mat;
    end
    
    % 6. Assign Right-Hand Side Vector
    beq(row) = b_val;
end