function [s,s_ter,u,du,y,err,...
    fi_s_min_x0,fi_s_max_x0,fi_s_ter_min_x0,fi_s_ter_max_x0,...
    fi_u_min_x0,fi_u_max_x0,fi_du_min_x0,fi_du_max_x0,...
    fi_y_min_x0,fi_y_max_x0,feas] = ...
    compute_box_constraints(x,mpc,check_feas)

%states
[s,s_ter] = get_x(x,mpc.nx,mpc.nu,mpc.N,mpc.Nx);
% control actions
u = get_u(x,mpc.nx,mpc.nu,mpc.N,mpc.Nu);
% differential control action
du = diff_u(u,u_prev,mpc.nu,mpc.N,mpc.Nu);
% system outputs
y = get_y(s,u,d,mpc.nx,mpc.nu,mpc.ny,mpc.nd,mpc.N,mpc.Ny,...
    mpc.C,mpc.D,mpc.Dd,mpc.Nd);
% error signal
err = get_error(r,y,mpc.N,mpc.Ny);

feas = 1;
% State box constraints
if ~isempty(mpc.x_min) && ~isempty(mpc.x_max)
    [fi_s_min_x0,fi_s_max_x0] = fi_box_fun(s,mpc.x_min,mpc.x_max,mpc.Nx,mpc.nx);
    if check_feas && ( any(fi_s_min_x0>0) || any(fi_s_max_x0>0) )
        feas = 0;
    end
end

% Terminal State box constraints
if ~isempty(mpc.x_ter_min) && ~isempty(mpc.x_ter_max) && feas
    [fi_s_ter_min_x0,fi_s_ter_max_x0] = fi_box_fun(s_ter,mpc.x_ter_min,mpc.x_ter_max,mpc.nx,mpc.nx);
    if check_feas && ( any(fi_s_ter_min_x0>0) || any(fi_s_ter_max_x0>0) )
        feas = 0;
    end
end

% Control box constraints
if ~isempty(mpc.u_min) && ~isempty(mpc.u_max) && feas
    [fi_u_min_x0,fi_u_max_x0] = fi_box_fun(u,mpc.u_min,mpc.u_max,mpc.Nu,mpc.nu);
    if check_feas && ( any(fi_u_min_x0>0) || any(fi_u_max_x0>0) )
        feas = 0;
    end
end

% Differential Control box constraints
if ~isempty(mpc.du_min) && ~isempty(mpc.du_max) && feas
    [fi_du_min_x0,fi_du_max_x0] = fi_box_fun(du,mpc.du_min,mpc.du_max,mpc.Nu,mpc.nu);
    if check_feas && ( any(fi_du_min_x0>0) || any(fi_du_max_x0>0) )
        feas = 0;
    end
end

% Outputs box constraints
if ~isempty(mpc.y_min) && ~isempty(mpc.y_max) && feas
    [fi_y_min_x0,fi_y_max_x0] = fi_box_fun(y,mpc.y_min,mpc.y_max,mpc.Ny,mpc.ny);
    if check_feas && ( any(fi_y_min_x0>0) || any(fi_y_max_x0>0) )
        feas = 0;
    end
end

end