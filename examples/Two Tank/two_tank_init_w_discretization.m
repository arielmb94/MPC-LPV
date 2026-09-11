%% Parameters

% Tank area and gravity value
Ab = 1;
g = 9.81;

% Initial condition
h1 = 0.45;
h2 = 0.45;

% Sampling time
Ts = 0.01;

%% Create MPC object

N = 10;         % Prediction Horizon

% Create mpc struct
mpc = init_mpc(N);
%% LTI system

% Get LPV model frozen at current state vector (Continuous-Time)
Ac = [-sqrt(2*g)*sqrt(h1)/(Ab*h1) 0;
       sqrt(2*g)*sqrt(h1)/(Ab*h1) -sqrt(2*g)*sqrt(h2)/(Ab*h2)];
Bc = [1/Ab; 0];
Bdc = []; % No disturbance matrix in this example

% Discretize the system
% You can easily swap 'forward' for 'backward' or 'tustin'
[Ad, Bd, Bdd] = init_discretize_system(Ac, Bc, Bdc, Ts, 'tustin');

% Initialize system dynamics
mpc = init_mpc_dynamics(mpc, Ad, Bd, Bdd);

% Tracking objective is the water height on the second tank
C = [0 1];

mpc = init_mpc_output(mpc,C,[],[]);

%% Constraints

% State constraints
x_min = 0.01*ones(mpc.nx,1);
x_max = 1*ones(mpc.nx,1);
mpc = init_mpc_state_cnstr(mpc,x_min,x_max);

% Control input constraints
u_min = 0*ones(mpc.nu,1);
u_max = 10*ones(mpc.nu,1);
mpc = init_mpc_u_cnstr(mpc,u_min,u_max);

% Control inputs variation constraints
du_min = -0.1*ones(mpc.nu,1);
du_max = 0.1*ones(mpc.nu,1);
mpc = init_mpc_delta_u_cnstr(mpc,du_min,du_max);

% Output constraints
y_min = [];
y_max = [];
%mpc = init_mpc_output_cnstr(mpc,y_min,y_max);

%% General Linear Inequalities

% General Linear Inequalities are defined as:
% yh = Ch*x+Dh*u+Ddh*di
Ch = [];
Dh = [];
Ddh = [];

h_min = [];
h_max = [];

%mpc = init_mpc_lin_custom_cnstr(mpc,Ch,Dh,Ddh,h_min,h_max);

%% Terminal Ingredients

% Terminal ingredients are computed using the dLQR method
Qx = diag([30 30]);         % State Penalty
Ru = 1;                     % Control Penalty
x_ref_is_y = 0;             % The terminal reference cannnot be extracted 
                            % from the mpc tracking reference

% Initialize terminal ingredients using the dLQR method                           
[mpc] = init_mpc_ter_ingredients_dlqr(mpc,Qx,Ru,x_ref_is_y);

%% Costs

% Tracking penalty
Qe = diag(50*ones(mpc.ny,1));
mpc = init_mpc_Tracking_cost(mpc,Qe);

% Control inputs variation penalty
Rdu = 1;
mpc = init_mpc_ControlRate_cost(mpc,Rdu);

% Control penalty
Ru = [];    % Quadratic penalty on control action u'*Ru*u
ru = [];    % Linear penalty on control action vector: ru'*u 
%mpc = init_mpc_Control_cost(mpc,Ru,ru);

%% Performance Cost Matrix

% Performance Vector are defined as:
% z = Cz*x+Dz*u+Ddz*dz
Cz = [];
Dz = [];
Dsuz = [];
Ddz = [];

Qz = [];    % Quadratic penalty on performance vector: z'*Qz*z
qz = [];    % Linear penalty on performance vector: qz'*z 

% Init performance cost
%mpc = init_mpc_Lin_Custom_cost(mpc,Cz,Dz,Dsuz,Ddz,Qz,qz);

%% Init conditions for simulation

% use warm start function to get optimization vector initial value
x_prev = [h1; h2];
u_prev = 3.7;
mpc = build_chronos_mpc(mpc,x_prev,u_prev);