%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Course: AAE43900 - 001

% Date: November 9, 2024

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

wind_velocity = 4; % [m/s]
rail_length = 1;   % [m]
C_D = 0.629; % [1] drag coefficient
rho = 1.225; % [kg/m^3] air density
g = 9.81; % [m/s^2]
A_rel = pi * (0.07874/2)^2; % [m^2] cross-sectional area
tb = 1.2; % [s] burnout time, pulled from given plot of thrust vs. time:
% https://wildmanrocketry.com/products/g74-6w?_pos=3&_sid=99a5eab6f&_ss=r
t_chute = tb + 6.5; % [s] time chute fully deploys after burnout
chute_CdA = 0.16; % [m^2] product of chute drag coefficient and area

constants = [wind_velocity; rail_length; C_D; rho; g; A_rel; tb; ...
    t_chute; chute_CdA];

% time constants
tstep = 0.01; % [s]
t_f = 40;    % [s]

% initial values
velocity_0 = 0; % [m/s]
theta_0 = 83.42 * pi/180;    % [rad]
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
subplot(2, 1, 1)
plot(z, h); grid on;
xlabel("Range [m]")
ylabel("Altitude [m]")
xlim([0 600])
ylim([0 200])

subplot(2, 1, 2)
plot(t, theta); grid on;
xlabel("Time [s]")
ylabel("Theta [rad]")
set(gcf, 'Position', [500 100 700 800])

figure(2);
plot(t, z);
title("range")

figure(3);
plot(t, h);
title("altitude")


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
    C_D = constants(3); % [1] drag coefficient
    rho = constants(4); % [kg/m^3] density of air
    g = constants(5);   % [m/s^2] constant of gravitation
    A_rel = constants(6); % [m^2] area of cross-section
    tb = constants(7); % [s] burnout time
    t_chute = constants(8); % [s] chute deployment time
    chute_CdA = constants(9); % [m^2] chute Cd * A
    
    % f(v, theta)
    V_rel = sqrt((velocity * sin(theta))^2 + (wind_velocity + ...
        velocity * cos(theta))^2); % [m/s] relative velocity
    psi = atan(velocity * sin(theta) / (wind_velocity + velocity * ...
        cos(theta))); % [rad] flight path angle
    % psi = acos((wind_velocity + velocity*cos(theta)) / V_rel);
    % psi = asin((velocity*sin(theta) / V_rel));
    
    % f(V_rel)
    if(t > t_chute)
        % [N] drag force with chute
        D = 0.5 * rho * V_rel^2 * (C_D * A_rel + chute_CdA);
    else
        D = 0.5 * rho * V_rel^2 * C_D * A_rel; % [N] drag force, no chute
    end

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
    trajectory_fun(1) = ((F - D) * cos(psi - theta) - ...
        m * g * sin(theta)) / m;
    
    dist = sqrt(altitude^2 + range^2);    
    % [rad/s] dtheta/dt
    if(dist < rail_length)
        trajectory_fun(2) = 0;
    else
        trajectory_fun(2) = ((F - D) * sin(psi - theta) - ...
            m * g * cos(theta)) / (m * velocity);
    end
    
    % [m/s] dh/dt
    trajectory_fun(3) = velocity * sin(theta);
    
    % [m/s] dz/dt
    trajectory_fun(4) = velocity * cos(theta);

end

% [N] thrust function
function thrust = thrust(time)
    thrust = 78; % [N]
end

% [kg] mass function
function mass = mass(time)
    mass = 0.921 - 0.087 * time; % [kg]
end