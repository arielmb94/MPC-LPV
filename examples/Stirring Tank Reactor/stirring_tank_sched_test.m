%% Call the mpc problem initialization script
stirring_tank_init

%% Define simulation duration and reference parameters
tsim = 240; % seconds
Sim_samples = tsim/Ts;
time = 0:Ts:tsim-Ts;

% Define concentration state reference
clear ref_c_vec
ref_c_vec = zeros(1,Sim_samples);
ref_c_vec(time < 90) = 0.27 + (0.65 - 0.27) * time(time < 90) / 90;
ref_c_vec(time >= 90 & time < 180) = 0.65;
ref_c_vec(time >= 180) = 0.65 - (0.65 - 0.3) * (time(time >= 180) - 180) / 60;
xref = [ref_c_vec(1); -M/(log(1/(theta_f*k*ref_c_vec(1))*(1-ref_c_vec(1))))];

%% Setup Functions and Bounds
my_sched_fun = @(x) calc_rho_stirring_tank(x, k, M);
my_jacob_fun = @(x) calc_jacob_stirring_tank(x, k, M);

% Physical bounds for recursive clipping
x_min = [0.0; 0.1];
x_max = [1.0; 100.0];
rho_min = my_sched_fun(x_min);
rho_max = my_sched_fun(x_max);

% Clear storage variables
clear c_dat v_dat u_dat ref_v_dat t_dat
n_rho = 3;
n_iter = 5; 

% --- DEFINE TARGET INSTANT FOR ANALYSIS ---
target_k = 5;

% Vectors to store predictions at the target instant
Pk_frozen_target = zeros(mpc.N * n_rho, 1);
Pk_iter_target = zeros(mpc.N * n_rho, 1);
Pk_iter_fast_target = zeros(mpc.N * n_rho, 1);
Pk_rec_target = zeros(mpc.N * n_rho, 1);

%% Run Simulation Loop
for i = 1:Sim_samples
    ref_c = ref_c_vec(i);
    ref_v = -M/(log(1/(theta_f*k*ref_c)*(1-ref_c)));
    ck = x_prev(1);
    vk = x_prev(2);
    
    tic
    
    % LPV Affine Matrices
    A0 = eye(2) + Ts * [-1/theta_f, 0; 0, -1/theta_f];
    A1 = Ts * [-1, 0; 1, 0];
    A2 = Ts * [ 0, -1; 0, 0];
    A3 = zeros(2,2);
    B0 = Ts * [0; alpha*xc];
    B1 = zeros(2,1);
    B2 = zeros(2,1);
    B3 = Ts * [0; -alpha];
    Bd0 = Ts * [1/theta_f, 0; xf/theta_f, 0];
    Bd1 = zeros(2,2);
    Bd2 = Ts * [0, 1; 0, 0];
    Bd3 = zeros(2,2);
    d = [1;vk];
    xref = [ref_c;ref_v];
    
    % --- COMPARE 4 METHODS AT TARGET INSTANT ---
    if i == target_k
        Pk_frozen_target = compute_schedul_frozen(mpc, x_prev, my_sched_fun);
        
        [Pk_iter_target, ~] = compute_schedul_iterative(mpc, x0, x_prev, u_prev, xref, d, [], ...
            my_sched_fun, n_rho, n_iter, ...
            A0, B0, Bd0, A1, A2, A3, B1, B2, B3, Bd1, Bd2, Bd3);
            
        Pk_iter_fast_target = compute_schedul_iterative_fast(mpc, x0, x_prev, my_sched_fun, n_rho);
        
        Pk_rec_target = compute_schedul_recursive(mpc, x0, x_prev, n_rho, ...
            my_sched_fun, my_jacob_fun, rho_min, rho_max);
    end
    
    % --- CONTINUE SIMULATION (Defaulting to Recursive Method) ---
    Pk = compute_schedul_recursive(mpc, x0, x_prev, n_rho, my_sched_fun, my_jacob_fun, rho_min, rho_max);
    
    mpc = update_mpc_sys_dynamics(mpc, A0, B0, Bd0, Pk, n_rho, A1, A2, A3, B1, B2, B3, Bd1, Bd2, Bd3);
    
    [u_prev, x0] = mpc_solve(mpc, x0, x_prev, u_prev, xref, d, [], [], [], ...
        A0, Bd0, Pk, n_rho, A1, A2, A3, B1, B2, B3, Bd1, Bd2, Bd3);
        
    tk = toc;
    
    % Store variables
    c_dat(:,i) = ck;
    v_dat(:,i) = vk;
    u_dat(i) = u_prev;
    ref_v_dat(i) = ref_v;
    t_dat(i) = tk;
    
    % Forward Euler step of Stirring Tank nonlinear dynamics
    ck = ck + Ts*((1-ck)/theta_f - k*ck*exp(-M/vk));
    vk = vk + Ts*((xf-vk)/theta_f + k*ck*exp(-M/vk)-alpha*u_prev*(vk-xc));
    x_prev = [ck;vk];
    
    % Shift x0 vector for warm-starting the next iteration
    % x0 = shift_warm_start(x0, mpc);
end

%% Post-Processing: Extract Real Pk for Target Instant
Pk_real_target = zeros(mpc.N * n_rho, 1);
for j = 1:mpc.N
    idx = target_k + j - 1;
    if idx > Sim_samples
        idx = Sim_samples; 
    end
    x_real = [c_dat(idx); v_dat(idx)];
    Pk_real_target((j-1)*n_rho + 1 : j*n_rho) = my_sched_fun(x_real);
end

%% Extract Individual Scheduling Variables for Plotting
% N x 5 Matrices (Col 1=Real, 2=Frozen, 3=Iterative, 4=Iterative Fast, 5=Recursive)
rho1_plot = [Pk_real_target(1:n_rho:end), Pk_frozen_target(1:n_rho:end), Pk_iter_target(1:n_rho:end), Pk_iter_fast_target(1:n_rho:end), Pk_rec_target(1:n_rho:end)];
rho2_plot = [Pk_real_target(2:n_rho:end), Pk_frozen_target(2:n_rho:end), Pk_iter_target(2:n_rho:end), Pk_iter_fast_target(2:n_rho:end), Pk_rec_target(2:n_rho:end)];
rho3_plot = [Pk_real_target(3:n_rho:end), Pk_frozen_target(3:n_rho:end), Pk_iter_target(3:n_rho:end), Pk_iter_fast_target(3:n_rho:end), Pk_rec_target(3:n_rho:end)];
k_steps = 1:mpc.N;

%% Target Instant Comparative Plot
figure('Name', ['Scheduling Parameters Estimation at k=' num2str(target_k)], 'Color', 'w', 'Position', [100, 100, 800, 700]);

% --- Subplot 1: Rho 1 ---
subplot(3,1,1);
plot(k_steps, rho1_plot(:,1), 'k-', 'LineWidth', 2.5); hold on;
plot(k_steps, rho1_plot(:,2), '-', 'Color', [0.7 0.7 0.7], 'LineWidth', 1.5); 
plot(k_steps, rho1_plot(:,3), 'r-.', 'LineWidth', 1.5); 
plot(k_steps, rho1_plot(:,4), 'm:', 'LineWidth', 2); 
plot(k_steps, rho1_plot(:,5), 'b--', 'LineWidth', 1.5); 
grid on;
ylabel('\rho_1');
title(['Scheduling parameters and corresponding estimates at k=' num2str(target_k)]);
legend('Real Trajectory', 'Frozen', 'Iterative SQP', 'Iterative Fast', 'Recursive', 'Location', 'best');

% --- Subplot 2: Rho 2 ---
subplot(3,1,2);
plot(k_steps, rho2_plot(:,1), 'k-', 'LineWidth', 2.5); hold on;
plot(k_steps, rho2_plot(:,2), '-', 'Color', [0.7 0.7 0.7], 'LineWidth', 1.5);
plot(k_steps, rho2_plot(:,3), 'r-.', 'LineWidth', 1.5);
plot(k_steps, rho2_plot(:,4), 'm:', 'LineWidth', 2);
plot(k_steps, rho2_plot(:,5), 'b--', 'LineWidth', 1.5);
grid on;
ylabel('\rho_2');

% --- Subplot 3: Rho 3 ---
subplot(3,1,3);
plot(k_steps, rho3_plot(:,1), 'k-', 'LineWidth', 2.5); hold on;
plot(k_steps, rho3_plot(:,2), '-', 'Color', [0.7 0.7 0.7], 'LineWidth', 1.5);
plot(k_steps, rho3_plot(:,3), 'r-.', 'LineWidth', 1.5);
plot(k_steps, rho3_plot(:,4), 'm:', 'LineWidth', 2);
plot(k_steps, rho3_plot(:,5), 'b--', 'LineWidth', 1.5);
grid on;
ylabel('\rho_3');
xlabel('Prediction Horizon Step (j)');