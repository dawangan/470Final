close all;
clc;
clear all;
% Define symbolic variables for angular velocities and Euler angles
syms wx wy wz real;
syms phi theta psi real;
syms Ix Iy Iz real;
%% translational variables
syms px py pz real;
syms vx vy vz real;
% External forces and torques as symbolic variables
tau_ext = sym('tau_ext', [3, 1], 'real');
force_ext = sym('force_ext', [3, 1], 'real');
%% define input u
u = [force_ext; tau_ext];
% Define the transformation matrix from angular velocities to Euler angle rates
T = [1, sin(phi)*tan(theta), cos(phi)*tan(theta);
     0, cos(phi), -sin(phi);
     0, sin(phi)/cos(theta), cos(phi)/cos(theta)];
T_phi = T(1,1:3);
T_theta = T(2,1:3);
T_psi = T(3,1:3);
% Define the symbolic inertia matrix I and angular velocity vector omega
I = diag([Ix, Iy, Iz]);
omega = [wx; wy; wz];
%% Translational motion dynamics definition
radius = 7071e3; % in z direction 
M_earth = 5.972e24;
M_satellite = 8211;
G_constant = 6.67e-11;
r= sqrt(px^2 + py^2 + pz^2);
Force_in_x = -G_constant * M_satellite * M_earth / (r^3) * px + u(1);
Force_in_y = -G_constant * M_satellite * M_earth / (r^3) * py + u(2);
Force_in_z = -G_constant * M_satellite * M_earth / (r^3) * pz + u(3);
F_position_translational = [vx; vy; vz]; % v
F_velocity_translational = [Force_in_x/M_satellite; Force_in_y/M_satellite; Force_in_z/M_satellite]; % a
% Angular velocity dynamics 
F_omega = I \ (-cross(omega, I * omega) + tau_ext);
% Euler angles dynamics in terms of angular velocities
F_euler = [T_phi; T_theta; T_psi] * omega;
% Combine angular velocity and Euler angle dynamics
F_full = [F_position_translational; F_velocity_translational; F_euler; F_omega];
% State vector includes both angular velocities and Euler angles
%% add translation
state_vector = [px; py; pz; vx; vy; vz; phi; theta; psi; wx; wy; wz;];
% Compute the full A matrix for the combined dynamics
A_full = jacobian(F_full, state_vector);
B_full = jacobian(F_full, u);
% Display the full A matrix
disp('Full A Matrix:');
disp(A_full);
% Display the full B matrix
disp('Full B Matrix:');
disp(B_full);
%%
% Define the output vector y as a symbolic variable
y = sym('y', [3, 1], 'real'); % Change the size of y according to your actual output vector
% Define the relationship between y and the state vector x
% For this example, assuming y is directly related to the Euler angles theta, psi
y = [theta; psi;];
% Define the input vector u as the external torques
u = tau_ext(1:3); % Assuming u corresponds to the external torques tau_ext
% Compute the C matrix by taking the Jacobian of the output vector with respect to the state vector
C = jacobian(y, state_vector);
% Compute the D matrix by taking the Jacobian of the output vector with respect to the input vector
D = jacobian(y, u);
% Display the C matrix
disp('C Matrix:');
disp(C);
% Display the D matrix
disp('D Matrix:');
disp(D);
%%
% Numerical values for the inertia matrix I and angular velocity omega_eq
I_vals = [4500, 5000, 5500]; % Inertia values [Ix, Iy, Iz] in kg·m²
tanslational_position_op = [radius, 0, 0];
translational_v_op = [0, 7502.4, 0]; 
theta_op = [0, 0, 0]; %Euler Angles Vector
omega_op = [0, 0, 0.0011]; %Angular velocities
% Equilibrium condition is when the derivative of omega is zero with zero external torques
% For a diagonal inertia matrix, this is satisfied when omega is aligned with an inertia axis
% Hence, we can define equilibrium points directly
omega_eq1 = [1, 0, 0]; % Spinning about the x-axis
omega_eq2 = [0, 1, 0]; % Spinning about the y-axis
omega_eq3 = [0, 0, 1]; % Spinning about the z-axis
% Time span for the simulation
tspan = [0, 100]; % Time span for the simulation
% A_numeric and B_numeric are calculated
A_numeric = double(subs(A_full, [Ix, Iy, Iz, px, py, pz, vx, vy, vz, phi, theta, psi, wx, wy, wz], [I_vals, tanslational_position_op, translational_v_op, theta_op, omega_op]));
B_numeric = double(subs(B_full, [Ix, Iy, Iz], I_vals));
% Calculate the eigenvalues of the A matrix
eigenvalues = eig(A_numeric);
% Print out the eigenvalues
disp('Eigenvalues of the system:');
disp(eigenvalues);
% Plot the eigenvalues
figure;
plot(real(eigenvalues), imag(eigenvalues), 'x');
xlabel('Real Part');
ylabel('Imaginary Part');
title('Eigenvalues of the Linearized System');
grid on;
axis equal;
%%
% Evaluate C and D numerically at the equilibrium point
C_numeric = double(subs(C, [Ix, Iy, Iz, phi, theta, psi, wx, wy, wz], [I_vals, 0, 0, 0, omega_op]));
D_numeric = double(subs(D, [Ix, Iy, Iz], I_vals));
disp('C Numeric Matrix:');
disp(C_numeric);
disp('D Numeric Matrix:');
disp(D_numeric);
% Construct the controllability matrix
ControllabilityMatrix = ctrb(A_numeric, B_numeric);
% Check if the system is controllable
if rank(ControllabilityMatrix) == size(A_numeric, 1)
    disp('The system is controllable.');
else
    disp('The system is not controllable.');
end
% Construct the observability matrix
ObservabilityMatrix = obsv(A_numeric, C_numeric);
% Check if the system is observable
if rank(ObservabilityMatrix) == size(A_numeric, 1)
    disp('The system is observable.');
else
    disp('The system is not observable.');
end
%%
% Adjusted numerical solver calls for full state dynamics
opts = odeset('RelTol',1e-6, 'AbsTol',1e-9);
u = [0;0;0];
% Solve the nonlinear dynamics numerically with full state dynamics
nonlinear_ode = @(t, state) nonlinear_dynamics_full(t, state, I_vals, [0, 0, 0]);
[time_nl, state_nl] = ode45(nonlinear_ode, tspan, initial_conditions, opts);
% % Calculate output y for nonlinear dynamics
% y_nl = arrayfun(@(idx) output_function(state_nl(idx, :).'), 1:size(state_nl,1), 'UniformOutput', false);
% y_nl = cell2mat(y_nl.');
% Solve the linear dynamics numerically with the full A matrix
linear_ode = @(t, state_vector) linear_dynamics_full(t, state_vector, u, A_numeric, B_numeric);
[time_l, state_l] = ode45(linear_ode, tspan, initial_conditions, opts);
% Directly calculate outputs after solving ODEs, using the state matrices
y_nl = C_numeric * state_nl(:, :).';  % Assuming state_nl has each state as a column
y_l = C_numeric * state_l(:, :).';    % Same assumption for state_l
% Ensure transposing is correctly done to match time vector dimensions
y_nl = y_nl.';  % Transpose to match time vector orientation
y_l = y_l.';    % Transpose to match time vector orientation
% Plotting results will now include Euler angles as well and output y
figure;
subplot(3,1,1);
plot(time_nl, state_nl(:,1:3), 'DisplayName', 'Nonlinear Dynamics');
hold on
plot(time_l, state_l(:,1:3), '--', 'DisplayName', 'Linear Dynamics');
xlabel('Time (s)');
ylabel('Euler Angles (rad)');
title('Euler Angles Dynamics');
legend('\phi', '\theta', '\psi');
grid on;
subplot(3,1,2);
plot(time_nl, state_nl(:,4:6), 'DisplayName', 'Nonlinear Dynamics');
hold on;
plot(time_l, state_l(:,4:6), '--', 'DisplayName', 'Linear Dynamics');
xlabel('Time (s)');
ylabel('Angular Velocity (rad/s)');
title('Angular Velocity Dynamics');
legend('wx', 'wy', 'wz');
grid on;
subplot(3,1,3);
plot(time_nl, y_nl, 'DisplayName', 'Nonlinear Output Dynamics');
hold on;
plot(time_l, y_l, '--', 'DisplayName', 'Linear Output Dynamics');
xlabel('Time (s)');
ylabel('Output (rad)');
title('Output Dynamics');
legend('\theta_{nl}', '\psi_{nl}','\theta_{l}', '\psi_{l}');
grid on;
ylim([-0.1,10])
% Update the dynamics functions to handle full state vector
function dstate_dt = linear_dynamics_full(t, state, u, A, B)
    dstate_dt = A * state + B * u; % state includes both Euler angles and angular velocities
end
function dstate_dt = nonlinear_dynamics_full(t, state, I_vals, tau_vals)
    % Extract angular velocities and Euler angles from the state vector
    phi = state(1);
    theta = state(2);
    psi = state(3);
    omega = state(4:6);
    
    % Compute the cross product term for angular velocity
    omega_cross = [0, -omega(3), omega(2); omega(3), 0, -omega(1); -omega(2), omega(1), 0];
    I = diag(I_vals);
    
    % Compute the dynamics for angular velocity
    domega_dt = I \ (- omega_cross * (I * omega) + tau_vals');
    
    % Transformation matrix from angular velocities to Euler angle rates
    T = [1, sin(phi)*tan(theta), cos(phi)*tan(theta);
         0, cos(phi), -sin(phi);
         0, sin(phi)/cos(theta), cos(phi)/cos(theta)];
    
    % Calculate Euler angles rates from angular velocities
    deuler_dt = T * omega;
    
    % Combine dynamics for Euler angles and angular velocities
    dstate_dt = [deuler_dt; domega_dt];
end