function Aeq = appendDynamics(Aeq, row, s_col, u_col, s_next_col, A_mat, B_mat, eye_next_mat)
% APPENDDYNAMICS Inserts a single block of system dynamics into Aeq.
% Optimized for MATLAB Coder (C Code Generation).

    % State mapping (A * s_k)
    if ~isempty(s_col) && ~isempty(A_mat)
        Aeq(row, s_col) = A_mat;
    end
    
    % Input mapping (B * u_k)
    if ~isempty(u_col) && ~isempty(B_mat)
        Aeq(row, u_col) = B_mat;
    end
    
    % Next state mapping (-I * s_{k+1})
    if ~isempty(s_next_col) && ~isempty(eye_next_mat)
        Aeq(row, s_next_col) = eye_next_mat;
    end
end