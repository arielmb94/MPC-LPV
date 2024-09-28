clear yk rk tk
r = -2*ones(ny,1);
tic
for k = 1:100

dist = randn*0.1;

if k>50
    r=2.5*ones(ny,1);
end

y = C*x_prev + D*u_prev + Dd*dist;


[u_prev,J,x0,t] = mpc_solve(x0,x_prev,u_prev,r,d,mpc,1e-2,1);
mpc.t = t;

x_prev = A*x_prev + B*u_prev + Bd*dist;


yk(:,k) = y;
rk(:,k) = r;
tk(k) = t;
Jk(k) = J;

end
toc/100
1/ans
J

close all
for y = 1:ny
    figure
    plot(1:k,yk(y,:),1:k,rk(y,:))
    grid on
end

figure
plot(1:k,tk)

%%

%states
s = get_x(x0,mpc.nx,mpc.nu,mpc.N,mpc.Nx)
% control actions
u = get_u(x0,mpc.nx,mpc.nu,mpc.N,mpc.Nu)
u0 = u(1:mpc.nu)
% differential control action
du = diff_u(u,u_prev,mpc.nu,mpc.N,mpc.Nu)
% system outputs
y = get_y(s,u,d,mpc.nx,mpc.nu,mpc.ny,mpc.nd,mpc.N,mpc.Ny,...
    mpc.C,mpc.D,mpc.Dd)
% error signal
err = get_error(r,y,mpc.N,mpc.Ny)

feas = 1;
% State box constraints
if ~isempty(mpc.x_min) & ~isempty(mpc.x_max)
    [fi_s_min_x0,fi_s_max_x0] = fi_box_fun(s,mpc.x_min,mpc.x_max,mpc.Nx,mpc.nx);
    if max(fi_s_min_x0)>0 || max(fi_s_max_x0)>0
        feas = 0;
    end
end

% Control box constraints
if ~isempty(mpc.u_min) & ~isempty(mpc.u_max) & feas
    [fi_u_min_x0,fi_u_max_x0] = fi_box_fun(u,mpc.u_min,mpc.u_max,mpc.Nu,mpc.nu);
    if max(fi_u_min_x0)>0 || max(fi_u_max_x0)>0
        feas = 0;
    end
end

% Differential Control box constraints
if ~isempty(mpc.du_min) & ~isempty(mpc.du_max) & feas
    [fi_du_min_x0,fi_du_max_x0] = fi_box_fun(du,mpc.du_min,mpc.du_max,mpc.Nu,mpc.nu);
    if max(fi_du_min_x0)>0 || max(fi_du_max_x0)>0
        feas = 0;
    end
end

% Outputs box constraints
if ~isempty(mpc.y_min) & ~isempty(mpc.y_max) & feas
    [fi_y_min_x0,fi_y_max_x0] = fi_box_fun(y,mpc.y_min,mpc.y_max,mpc.Ny,mpc.ny);
    if max(fi_y_min_x0)>0 || max(fi_y_max_x0)>0
        feas = 0;
    end
end

feas



