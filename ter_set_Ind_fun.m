function [fi_ter_x0,grad_ter,grad_ter_Ind_x0,hess_ter_Ind_x0] = ter_set_Ind_fun(s_ter,P,Nx,Nu,nx,nu,N)

% Terminal Set Constraint Evaluation
fi_ter_x0 = s_ter'*P*s_ter-1;

% Gradient of Terminal Set Constraint
grad_ter = P*s_ter;
% Gradient of Indicator Function for Terminal Set Constraint
grad_ter_Ind_x0 = zeros(Nx+Nu,1);
grad_ter_Ind_x0(nx*(N) + nu*(N+1) + 1 : end) = -grad_ter/fi_ter_x0;

% Hessian of Indicator Function for Terminal Set Constraint
hess_ter_Ind_x0 = zeros(Nx+Nu,Nx+Nu);
hess_ter_Ind_x0(nx*(N) + nu*(N+1) + 1 : end,nx*(N) + nu*(N+1) + 1 : end) = ...
    (grad_ter*grad_ter')/(fi_ter_x0^2) - P/fi_ter_x0;
 
end