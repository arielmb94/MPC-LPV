function [Ai, bi] = appendGeneralizedConstraint(Ai, bi, row, s_col, su_col, u_col, C_mat, Dsu_mat, D_mat, b_val)
% APPENDGENERALIZEDCONSTRAINT Inserts a unified linear constraint block.
% Form: C*s + D*u + Ddu*su + I*g - I*v = b_val

    % 1. State term (C * s_k)
    if ~isempty(s_col) && ~isempty(C_mat)
        Ai(row, s_col) = C_mat;
    end

    % 2. Rate/Past Input term (Dsu * su_k)
    if ~isempty(su_col) && ~isempty(Dsu_mat)
        Ai(row, su_col) = Dsu_mat;
    end
    
    % 3. Control Input term (D * u_k)
    if ~isempty(u_col) && ~isempty(D_mat)
        Ai(row, u_col) = D_mat;
    end
    
    % 4. Assign Right-Hand Side Vector
    bi(row) = b_val;
end