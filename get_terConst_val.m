function fi_ter_x0 = get_terConst_val(x_ref,s_ter,P)

% compute error at terminal state
x_err = x_ref-s_ter;

% Terminal Set Constraint Evaluation
fi_ter_x0 = x_err'*P*x_err-1;


end