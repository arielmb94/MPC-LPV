function t = init_t(x0,u_prev,C,D,Dd,d,sigma,x_low,x_up,u_low,u_up,du_low,du_up,y_low,y_up, ...
    N,Nx,Nu,Ny,nx,nu,ny,nd)

    s = get_x(x0,nx,nu,N,Nx);

    u = get_u(x0,nx,nu,N,Nu);

    du = diff_u(u,u_prev,nu,N,Nu);

    y = get_y(s,u,d,nx,nu,ny,nd,N,Ny,C,D,Dd);

    fi_sum = 0;
    ni = 0;

    if ~isempty(x_low) & ~isempty(x_up)

        [fi_s_low_x0,fi_s_up_x0] = fi_box_fun(s,x_low,x_up,Nx,nx);

        fi_sum = fi_sum + sum(-fi_s_low_x0) + sum(-fi_s_up_x0);
        ni = ni + length(fi_s_low_x0) + length(fi_s_up_x0);
    end

    if ~isempty(u_low) & ~isempty(u_up)
        
        [fi_u_low_x0,fi_u_up_x0] = fi_box_fun(u,u_low,u_up,Nu,nu);

        fi_sum = fi_sum + sum(-fi_u_low_x0) + sum(-fi_u_up_x0);
        ni = ni + length(fi_u_low_x0) + length(fi_u_up_x0);
    end

    if ~isempty(du_low) & ~isempty(du_up)
        
        [fi_du_low_x0,fi_du_up_x0] = fi_box_fun(du,du_low,du_up,Nu,nu);

        fi_sum = fi_sum + sum(-fi_du_low_x0) + sum(-fi_du_up_x0);
        ni = ni + length(fi_du_low_x0) + length(fi_du_up_x0);
    end

    if ~isempty(y_low) & ~isempty(y_up)
        
        [fi_y_low_x0,fi_y_up_x0] = fi_box_fun(y,y_low,y_up,Ny,ny);

        fi_sum = fi_sum + sum(-fi_y_low_x0) + sum(-fi_y_up_x0);
        ni = ni + length(fi_y_low_x0) + length(fi_y_up_x0);
    end
    
    t = sigma * fi_sum / ni;
    
end


 