% m is number of inequalities -> to be precomputed

function [u0,J,x] = mpc_solve(x0,x_prev,u_prev,r,d,mpc,eps,t)

    % number of variables
    n = length(x0);

    % number of equality constraints
    n_eq = size(mpc.Aeq,1); 

%     % Compute inequality functions at x0 to compute duality measure
%     t = init_t(x0,u_prev,mpc.C,mpc.D,mpc.Dd,d,sigma,mpc.x_min,mpc.x_max,...
%         mpc.u_min,mpc.u_max,mpc.du_min,mpc.du_max,mpc.y_min,mpc.y_max, ...
%         mpc.N,mpc.Nx,mpc.Nu,mpc.Ny,mpc.nx,mpc.nu,mpc.ny,mpc.nd);

    mpc.beq = update_beq(mpc.beq,mpc.A,x_prev,mpc.N,mpc.nx,mpc.Bd,d*0);

    x = x0;

    lambda2 = 1;
    while eps <= lambda2/2


        % Get variables

        %states
        s = get_x(x,mpc.nx,mpc.nu,mpc.N,mpc.Nx);
        % control actions
        u = get_u(x,mpc.nx,mpc.nu,mpc.N,mpc.Nu);
        % differential control action
        du = diff_u(u,u_prev,mpc.nu,mpc.N,mpc.Nu);
        % system outputs
        y = get_y(s,u,d,mpc.nx,mpc.nu,mpc.ny,mpc.nd,mpc.N,mpc.Ny,...
            mpc.C,mpc.D,mpc.Dd);
        % error signal
        err = get_error(r,y,mpc.N,mpc.Ny);

        % Compute gradient:

        % 1. Compute Values of each box Ineq. Function fi(x)

        % State box constraints
        if ~isempty(mpc.x_min) & ~isempty(mpc.x_max)
            [fi_s_min_x0,fi_s_max_x0] = fi_box_fun(s,mpc.x_min,mpc.x_max,mpc.Nx,mpc.nx);
        end

        % Control box constraints
        if ~isempty(mpc.u_min) & ~isempty(mpc.u_max)
            [fi_u_min_x0,fi_u_max_x0] = fi_box_fun(u,mpc.u_min,mpc.u_max,mpc.Nu,mpc.nu);
        end

        % Differential Control box constraints
        if ~isempty(mpc.du_min) & ~isempty(mpc.du_max)
            [fi_du_min_x0,fi_du_max_x0] = fi_box_fun(du,mpc.du_min,mpc.du_max,mpc.Nu,mpc.nu);
        end

        % Outputs box constraints
        if ~isempty(mpc.y_min) & ~isempty(mpc.y_max)
            [fi_y_min_x0,fi_y_max_x0] = fi_box_fun(y,mpc.y_min,mpc.y_max,mpc.Ny,mpc.ny);
        end


        % 2. Compute gradient of box inequalities at x0:

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

        % line search using backtracking
        %l = linesearch_backtracking(x0,@f0_fun,@fi_fun,grad_x0,delta_x_prim,t,m,alpha,beta);
        l = 1;
        % update x0
        x = x+l*delta_x_prim;
        %J = f0_fun_MPC(Qe,err,N,ny,R,du,nu,[],[]);

    end


    %states
    s = get_x(x,mpc.nx,mpc.nu,mpc.N,mpc.Nx);
    % control actions
    u = get_u(x,mpc.nx,mpc.nu,mpc.N,mpc.Nu);
    u0 = u(1:mpc.nu);
    % differential control action
    du = diff_u(u,u_prev,mpc.nu,mpc.N,mpc.Nu);
    % system outputs
    y = get_y(s,u,d,mpc.nx,mpc.nu,mpc.ny,mpc.nd,mpc.N,mpc.Ny,...
        mpc.C,mpc.D,mpc.Dd);
    % error signal
    err = get_error(r,y,mpc.N,mpc.Ny);

    J = f0_fun_MPC(mpc.Qe,err,mpc.N,mpc.ny,mpc.R,du,mpc.nu,[],[]);

end






