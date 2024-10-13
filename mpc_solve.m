% m is number of inequalities -> to be precomputed

function [u0,J,x,t] = mpc_solve(x0,x_prev,u_prev,r,d,mpc,eps)

    % number of variables
    n = length(x0);

    % number of equality constraints
    n_eq = size(mpc.Aeq,1); 

    mpc.beq = update_beq(mpc.beq,mpc.A,x_prev,mpc.N,mpc.nx,mpc.Bd,d,mpc.Nd);

    x = x0;
    forward_iter = 0;

    t = mpc.t;

    search_feas = false;

    while forward_iter<mpc.max_iter

        % for first iteration we assume x0 is feasible and wont check
        check_feas = false;
        % get mpc variables from optimimization vector x and constraint
        % information
        [s,s_ter,u,du,y,err,...
            fi_s_min_x0,fi_s_max_x0,fi_s_ter_min_x0,fi_s_ter_max_x0,...
            fi_u_min_x0,fi_u_max_x0,fi_du_min_x0,fi_du_max_x0,...
            fi_y_min_x0,fi_y_max_x0,feas] = ...
            get_state_constraint_info(x,u_prev,r,d,mpc,check_feas);
        % for following iteration check feasibility
        check_feas = true;

        lambda2 = 1;
        while eps <= lambda2/2


            % Compute gradient:

            % 1. Compute gradient of box inequalities at x0:
            % init inequalities gradient vector
            grad_fi_Ind = zeros(n,1);

            % state inequalities
            if ~isempty(mpc.x_min)
                grad_s_min_Ind_x0 = grad_box_Ind(fi_s_min_x0,mpc.gradXmin);

                grad_fi_Ind = grad_fi_Ind + grad_s_min_Ind_x0;
            end

            if ~isempty(mpc.x_max)
                grad_x_max_Ind_x0 = grad_box_Ind(fi_s_max_x0,mpc.gradXmax);

                grad_fi_Ind = grad_fi_Ind + grad_x_max_Ind_x0;
            end

            % terminal state inequalities
            if ~isempty(mpc.x_ter_min)
                grad_s_ter_min_Ind_x0 = grad_box_Ind(fi_s_ter_min_x0,mpc.gradXtermin);

                grad_fi_Ind = grad_fi_Ind + grad_s_ter_min_Ind_x0;
            end

            if ~isempty(mpc.x_ter_max)
                grad_x_ter_max_Ind_x0 = grad_box_Ind(fi_s_ter_max_x0,mpc.gradXtermax);

                grad_fi_Ind = grad_fi_Ind + grad_x_ter_max_Ind_x0;
            end

            % control inequalities
            if ~isempty(mpc.u_min)
                grad_u_min_Ind_x0 = grad_box_Ind(fi_u_min_x0,mpc.gradUmin);

                grad_fi_Ind = grad_fi_Ind + grad_u_min_Ind_x0;
            end

            if ~isempty(mpc.u_max)
                grad_u_max_Ind_x0 = grad_box_Ind(fi_u_max_x0,mpc.gradUmax);

                grad_fi_Ind = grad_fi_Ind + grad_u_max_Ind_x0;
            end

            % control differential inequalities
            if ~isempty(mpc.du_min)
                grad_du_min_Ind_x0 = grad_box_Ind(fi_du_min_x0,mpc.gradDeltaUmin);

                grad_fi_Ind = grad_fi_Ind + grad_du_min_Ind_x0;
            end

            if ~isempty(mpc.du_max)
                grad_du_max_Ind_x0 = grad_box_Ind(fi_du_max_x0,mpc.gradDeltaUmax);

                grad_fi_Ind = grad_fi_Ind + grad_du_max_Ind_x0;
            end

            % output inequalities
            if ~isempty(mpc.y_min)
                grad_y_min_Ind_x0 = grad_box_Ind(fi_y_min_x0,mpc.gradYmin);

                grad_fi_Ind = grad_fi_Ind + grad_y_min_Ind_x0;
            end

            if ~isempty(mpc.y_max)
                grad_y_max_Ind_x0 = grad_box_Ind(fi_y_max_x0,mpc.gradYmax);

                grad_fi_Ind = grad_fi_Ind + grad_y_max_Ind_x0;
            end

            % 4. Compute gradient of cost function at x0
            grad_f0 = grad_f0_MPC(mpc.Nx,mpc.Nu,mpc.gradErrQe,err,mpc.gradCtlrR,du,[],[],[]);
            % 5. Compute gradient at x0 : grad(J) = t*grad(f0)+grad(Phi)
            grad_J_x0 = t*grad_f0+grad_fi_Ind;


            % Compute Hessian Matrix

            % 1. Compute Hessian of box inequalities at x0:

            % init inequalities hessian vector
            hess_fi_Ind = zeros(n);

            % state inequalities
            if ~isempty(mpc.x_min)
                hess_s_min_Ind_x0 = hess_linear_Ind(fi_s_min_x0,mpc.hessXmin);

                hess_fi_Ind = hess_fi_Ind + hess_s_min_Ind_x0;
            end

            if ~isempty(mpc.x_max)
                hess_s_max_Ind_x0 = hess_linear_Ind(fi_s_max_x0,mpc.hessXmax);

                hess_fi_Ind = hess_fi_Ind + hess_s_max_Ind_x0;
            end

            % terminal state inequalities
            if ~isempty(mpc.x_ter_min)
                hess_s_ter_min_Ind_x0 = hess_linear_Ind(fi_s_ter_min_x0,mpc.hessXtermin);

                hess_fi_Ind = hess_fi_Ind + hess_s_ter_min_Ind_x0;
            end

            if ~isempty(mpc.x_ter_max)
                hess_s_ter_max_Ind_x0 = hess_linear_Ind(fi_s_ter_max_x0,mpc.hessXtermax);

                hess_fi_Ind = hess_fi_Ind + hess_s_ter_max_Ind_x0;
            end

            % control inequalities
            if ~isempty(mpc.u_min)
                hess_u_min_Ind_x0 = hess_linear_Ind(fi_u_min_x0,mpc.hessUmin);

                hess_fi_Ind = hess_fi_Ind + hess_u_min_Ind_x0;
            end

            if ~isempty(mpc.u_max)
                hess_u_max_Ind_x0 = hess_linear_Ind(fi_u_max_x0,mpc.hessUmax);

                hess_fi_Ind = hess_fi_Ind + hess_u_max_Ind_x0;
            end

            % control differential inequalities
            if ~isempty(mpc.du_min)
                hess_du_min_Ind_x0 = hess_linear_Ind(fi_du_min_x0,mpc.hessDeltaUmin);

                hess_fi_Ind = hess_fi_Ind + hess_du_min_Ind_x0;
            end

            if ~isempty(mpc.du_max)
                hess_du_max_Ind_x0 = hess_linear_Ind(fi_du_max_x0,mpc.hessDeltaUmax);

                hess_fi_Ind = hess_fi_Ind + hess_du_max_Ind_x0;
            end

            % output inequalities
            if ~isempty(mpc.y_min)
                hess_y_min_Ind_x0 = hess_linear_Ind(fi_y_min_x0,mpc.hessYmin);

                hess_fi_Ind = hess_fi_Ind + hess_y_min_Ind_x0;
            end

            if ~isempty(mpc.y_max)
                hess_y_max_Ind_x0 = hess_linear_Ind(fi_y_max_x0,mpc.hessYmax);

                hess_fi_Ind = hess_fi_Ind + hess_y_max_Ind_x0;
            end

            % 2. Compute Hessian of cost Function
            hess_f0 = mpc.hessCtrlTerm + mpc.hessErrTerm;
            % 3. Compute Hessian of f(x0,t):
            hess_J_x0 = t*hess_f0+hess_fi_Ind;

            % solve KKT system
            KKT = [hess_J_x0 mpc.Aeq';mpc.Aeq zeros(n_eq)];

            %delta_x = - linsolve(KKT,[grad_J_x0;Aeq*x-beq],opts);
            %delta_x = - linsolve(KKT,[grad_J_x0;zeros(n_eq,1)],opts);
            delta_x = - KKT\[grad_J_x0;mpc.Aeq*x-mpc.beq];
            delta_x_prim = delta_x(1:n);

            % compute lambda^2
            lambda2 = -grad_J_x0'*delta_x_prim;

            % Feasibility line search

            l = 1;
            xhat = x+l*delta_x_prim;

            [s,s_ter,u,du,y,err,...
            fi_s_min_x0,fi_s_max_x0,fi_s_ter_min_x0,fi_s_ter_max_x0,...
            fi_u_min_x0,fi_u_max_x0,fi_du_min_x0,fi_du_max_x0,...
            fi_y_min_x0,fi_y_max_x0,feas] = ...
            get_state_constraint_info(xhat,u_prev,r,d,mpc,check_feas);

            if feas
                x = xhat;
            else
                while ~feas
                    l = l*mpc.Beta;

                    xhat = x+l*delta_x_prim;

                    [s,s_ter,u,du,y,err,...
                        fi_s_min_x0,fi_s_max_x0,fi_s_ter_min_x0,fi_s_ter_max_x0,...
                        fi_u_min_x0,fi_u_max_x0,fi_du_min_x0,fi_du_max_x0,...
                        fi_y_min_x0,fi_y_max_x0,feas] = ...
                        get_state_constraint_info(xhat,u_prev,r,d,mpc,check_feas);
                end
                x = xhat;
                %t = t / mpc.eta_bck;
                search_feas = true;
            end

        end

%         check_feas = true;
% 
%         [s,s_ter,u,du,y,err,...
%             fi_s_min_x0,fi_s_max_x0,fi_s_ter_min_x0,fi_s_ter_max_x0,...
%             fi_u_min_x0,fi_u_max_x0,fi_du_min_x0,fi_du_max_x0,...
%             fi_y_min_x0,fi_y_max_x0,feas] = ...
%             get_state_constraint_info(xhat,u_prev,r,d,mpc,check_feas);

        u0 = u(1:mpc.nu);
                
        if search_feas
            forward_iter = mpc.max_iter;
            %t = t / mpc.eta_bck;
        else
            if t < mpc.t_max
                t = t * mpc.eta_fwd;
                if t > mpc.t_max
                    t = mpc.t_max;
                end
            end
            forward_iter = forward_iter + 1;
        end
        
    end

    J = f0_fun_MPC(mpc.Qe,err,mpc.N,mpc.ny,mpc.R,du,mpc.nu,[],[]);

end






