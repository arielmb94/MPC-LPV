%% Call the mpc problem initialization script
two_tank_init

%% Define simulation duration and reference parameters
tsim = 5; % seconds
Sim_samples = tsim/Ts;
time = 0:Ts:tsim-Ts;

% Define step reference
clear r
r = zeros(1,Sim_samples);
r(time<=2.5) = 0.7;
r(time>2.5) = 0.25;

% Low-pass reference filter
tau = 0.1;      
xf = h2;        

%% Setup Functions and Bounds
n_rho = 2;
n_iter = 20; 

% LPV Functions
my_sched_fun = @(x) [1/sqrt(max(x(1), 1e-4)); 
                     1/sqrt(max(x(2), 1e-4))];
                     
my_jacob_fun = @(x) [-0.5/(max(x(1), 1e-4)^1.5), 0; 
                      0, -0.5/(max(x(2), 1e-4)^1.5)];

% Physical bounds (Based on x_min = 0.01 and x_max = 1)
% Max and min are inverted because 1/sqrt(smaller_value) yields the LARGEST rho.
rho_min = my_sched_fun([1.0; 1.0]); 
rho_max = my_sched_fun([0.01; 0.01]); 

% Clear storage variables
clear rf_dat rk h1_dat h2_dat u_dat t_dat

% --- DEFINE THE ANALYSIS INSTANT HERE ---
%target_k = round(2.6 / Ts); % Instant right after the step drop
target_k = 5;

% Vectors to store predictions at the target instant
Pk_frozen_target    = zeros(mpc.N * n_rho, 1);
Pk_iter_target      = zeros(mpc.N * n_rho, 1);
Pk_iter_fast_target = zeros(mpc.N * n_rho, 1);
Pk_rec_target       = zeros(mpc.N * n_rho, 1);

%% Run Simulation Loop
for k = 1:Sim_samples   
    h1  = x_prev(1);        
    h2  = x_prev(2);        
    
    % Low pass reference filter step
    xf = xf + Ts*(-xf/tau+r(k)/tau);
    x_ref = [xf;xf];
    
    tic;
    
    % --- LPV AFFINE MATRICES FOR THE TWO TANK ---
    k_tank = sqrt(2*g)/Ab;
    A0 = eye(2);
    A1 = Ts * [-k_tank, 0; 
                k_tank, 0];
    A2 = Ts * [0, 0; 
               0, -k_tank];
               
    B0 = mpc.B; % mpc.B is already Ts * [1/Ab; 0]
    B1 = zeros(2,1);
    B2 = zeros(2,1);
    
    % --- COMPARISON TEST OF THE 4 METHODS AT THE TARGET INSTANT ---
    if k == target_k
        Pk_frozen_target = compute_schedul_frozen(mpc, x_prev, my_sched_fun);
        
        % New call with x_ref and d=[]
        [Pk_iter_target, ~] = compute_schedul_iterative(mpc, x0, x_prev, u_prev, xf, [], x_ref, ...
                                              my_sched_fun, n_rho, n_iter, ...
                                              A0, B0, [], A1, A2, B1, B2);
                                              
        Pk_iter_fast_target = compute_schedul_iterative_fast(mpc, x0, x_prev, my_sched_fun, n_rho);
                                              
        Pk_rec_target = compute_schedul_recursive(mpc, x0, x_prev, n_rho, ...
                                             my_sched_fun, my_jacob_fun, rho_min, rho_max);
    end
    
    % --- SIMULATION CONTINUES (Using Iterative Fast as default) ---
    Pk = compute_schedul_iterative_fast(mpc, x0, x_prev, my_sched_fun, n_rho);
    
    % Dynamic update WITHOUT Bd matrices
    mpc = update_mpc_sys_dynamics(mpc, A0, B0, [], Pk, n_rho, A1, A2, B1, B2);
    
    [u_prev,x0] = mpc_solve(mpc, x0, x_prev, u_prev, xf, [], x_ref, [], [], ...
                            A0, [], Pk, n_rho, A1, A2, B1, B2);
    tk = toc;
    
    % Store variables
    rf_dat(:,k) = xf;
    h1_dat(:,k) = h1;
    h2_dat(:,k) = h2;
    u_dat(k) = u_prev;
    t_dat(k) = tk;
    
    % Forward Euler step
    h1 = h1 + Ts*(u_prev/Ab-sqrt(2*g)*sqrt(h1)/Ab);
    h2 = h2 + Ts*(sqrt(2*g)*sqrt(h1)/Ab-sqrt(2*g)*sqrt(h2)/Ab);
    
    x_prev = [h1;h2];
    
    % Perform shift for the next iteration
    % x0 = shift_warm_start(x0, mpc);
end

%% Post-Processing and Plots
close all
Pk_real_target = zeros(mpc.N * n_rho, 1);
for j = 1:mpc.N
    idx = target_k + j - 1; 
    if idx > Sim_samples
        idx = Sim_samples; 
    end
    
    x_real = [h1_dat(idx); h2_dat(idx)];
    Pk_real_target((j-1)*n_rho + 1 : j*n_rho) = my_sched_fun(x_real);
end

rho1_plot = [Pk_real_target(1:n_rho:end), Pk_frozen_target(1:n_rho:end), Pk_iter_target(1:n_rho:end), Pk_iter_fast_target(1:n_rho:end), Pk_rec_target(1:n_rho:end)];
rho2_plot = [Pk_real_target(2:n_rho:end), Pk_frozen_target(2:n_rho:end), Pk_iter_target(2:n_rho:end), Pk_iter_fast_target(2:n_rho:end), Pk_rec_target(2:n_rho:end)];
k_steps = 1:mpc.N;

figure('Name', ['Two Tank Estimation at k=' num2str(target_k)], 'Color', 'w', 'Position', [100, 100, 800, 500]);
subplot(2,1,1);
plot(k_steps, rho1_plot(:,1), 'k-', 'LineWidth', 2.5); hold on;
plot(k_steps, rho1_plot(:,2), '-', 'Color', [0.7 0.7 0.7], 'LineWidth', 1.5);
plot(k_steps, rho1_plot(:,3), 'r-.', 'LineWidth', 1.5);
plot(k_steps, rho1_plot(:,4), 'm:', 'LineWidth', 2);
plot(k_steps, rho1_plot(:,5), 'b--', 'LineWidth', 1.5);
grid on; ylabel('\rho_1');
title(['Scheduling parameters estimates at k=' num2str(target_k) ' (Two Tank)']);
legend('Real Trajectory', 'Frozen', 'Iterative SQP', 'Iterative Fast', 'Recursive', 'Location', 'best');

subplot(2,1,2);
plot(k_steps, rho2_plot(:,1), 'k-', 'LineWidth', 2.5); hold on;
plot(k_steps, rho2_plot(:,2), '-', 'Color', [0.7 0.7 0.7], 'LineWidth', 1.5);
plot(k_steps, rho2_plot(:,3), 'r-.', 'LineWidth', 1.5);
plot(k_steps, rho2_plot(:,4), 'm:', 'LineWidth', 2);
plot(k_steps, rho2_plot(:,5), 'b--', 'LineWidth', 1.5);
grid on; ylabel('\rho_2'); xlabel('Prediction Horizon Step (j)');

%% Standard Control Plots
figure
ax1 = subplot(2,1,1);
plot(time,r,'r',time,rf_dat,'--r',time,h1_dat,'g',time,h2_dat,'b')
grid on; legend('Reference','Filtered Reference','h1','h2'); xlabel('Time (s)'); ylabel('Water Height');
ax2 = subplot(2,1,2);
plot(time,u_dat,time(1:end-1),diff(u_dat))
grid on; legend('u','\Delta u'); xlabel('Time (s)'); ylabel('Input Mass Flow');
linkaxes([ax1,ax2],'x');