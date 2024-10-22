%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Course: AAE43900 - 001

% Date: October 27, 2024

% Co-Authors:
% Daniel DeConti
% Ryan Morey

% Rocket Project

% Description: model the trajectory of a small rocket with a solid
% propellant motor.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
format shortg
format compact
close all
clear
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs/Constants:

wind_velocity = 6; % [m/s]
rail_length = 0.5;   % [m]
C_L = 1.5;
C_D = 0.01;
rho = 1.225; % [kg/m^3] air density
g = 9.81; % [m/s^2]
A_rel = pi * (0.07874/2)^2; % [m^2] cross-sectional area
tb = 1.2; % [s] burnout time, pulled from given plot of thrust vs. time:
% https://wildmanrocketry.com/products/g74-6w?_pos=3&_sid=99a5eab6f&_ss=r

constants = [wind_velocity; rail_length; C_L; C_D; rho; g; A_rel; tb];

% time constants
tstep = 1; % [s]
t_f = 53;    % [s]

% initial values
velocity_0 = 0; % [m/s]
theta_0 = pi/4;    % [rad]
altitude_0 = 0; % [m]
range_0 = 0;    % [m]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% call ode45 for trajectory_fun()

tr = 0:tstep:t_f; % specify timesteps [s]

op = odeset("RelTol", 1e-6, "AbsTol", 1e-6); % set reasonable tolerances

% numerically integrate using state variable vector for v(t), theta(t),
% h(t), and z(t) for given drag inputs
[t, trajectory_vars] = ode45(@(t, trajectory_vars) trajectory_fun(t, ...
    trajectory_vars, constants), tr, [velocity_0; theta_0; altitude_0; ...
    range_0], op);

v = trajectory_vars(:,1);
theta = trajectory_vars(:,2);
h = trajectory_vars(:,3);
z = trajectory_vars(:,4);
figure(1)
plot(z, h); grid on;
xlabel("Range [m]")
ylabel("Altitude [m]")

figure(2)
plot(t, theta)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ODE45 solver function for velocity, theta, altitude, and range
function trajectory_fun = trajectory_fun(t, trajectory_vars, constants)

    % unpack trajectory variables
    velocity = trajectory_vars(1); % [m/s] total velocity
    theta = trajectory_vars(2);    % [rad] angle of total velocity
    altitude = trajectory_vars(3); % [m] altitude in inertial frame
    range = trajectory_vars(4);    % [m] range in inertial frame
    
    % unpack constants
    wind_velocity = constants(1); % [m/s] assume constant wind
    rail_length = constants(2);   % [m/s] length of rail
    C_L = constants(3); % [1] lift coefficient
    C_D = constants(4); % [1] drag coefficient
    rho = constants(5); % [kg/m^3] density of air
    g = constants(6);   % [m/s^2] constant of gravitation
    A_rel = constants(7); % [m^2] area of cross-section
    tb = constants(8); % [s] burnout time
    
    % f(v, theta)
    psi = atan(velocity * sin(theta) / (wind_velocity + velocity * ...
        cos(theta))); % [rad] flight path angle
    V_rel = sqrt((velocity * sin(theta))^2 + (wind_velocity + ...
        velocity * cos(theta))^2); % [m/s] relative velocity
    
    % f(V_rel)
    L = 0.5 * rho * V_rel^2 * C_L * A_rel; % [N] lift force
    D = 0.5 * rho * V_rel^2 * C_D * A_rel; % [N] drag force

    % f(t)
    if(t < tb) % before burnout
        m = mass(t); % [kg] mass of rocket at current time
        F = thrust(t); % [N] thrust force
    else % after burnout
        m = 0.921 - 0.087 * tb; % [kg] mass of rocket
        F = 0; % [N] thrust force
    end
    
    % return values
    trajectory_fun = linspace(1, 1, 4)';
    % [m/s^2] dv/dt
    trajectory_fun(1) = ((F - D) * cos(psi - theta) - L * sin(psi - ...
        theta) - m * g * sin(theta)) / m;
    
    dist = sqrt(altitude^2 + range^2);
    % [rad/s] dtheta/dt
    if(dist < rail_length)
        trajectory_fun(2) = 0;
    else
        trajectory_fun(2) = ((F - D) * sin(psi - theta) + L * cos(psi - ...
            theta) - m * g * cos(theta)) / (m * velocity);
    end
    
    % [m/s] dh/dt
    trajectory_fun(3) = velocity * sin(theta);
    % trajectory_fun(3) = trajectory_fun(1) * sin(trajectory_fun(2));
    
    % [m/s] dz/dt
    trajectory_fun(4) = velocity * cos(theta);
    % trajectory_fun(4) = trajectory_fun(1) * cos(trajectory_fun(2));
end

% [N] thrust function
function thrust = thrust(time)
    thrust = 78; % [N]
end

% [kg] mass function
function mass = mass(time)
    mass = 0.921 - 0.087 * time; % [kg]
end