%#ok<*NBRAK2>

% Commented out tests for other program

    % if(dist < rail_length)
    %     V_rel = sqrt((velocity * sin(theta))^2 + (0 + ...
    %         velocity * cos(theta))^2); % [m/s] relative velocity
    % else
    %     V_rel = sqrt((velocity * sin(theta))^2 + (wind_velocity + ...
    %         velocity * cos(theta))^2); % [m/s] relative velocity
    % end

%%
close all;clear;clc;

Vw = 6; % [m/s]
Lrod = 0.5;   % [m]
C_L = 1.5;
C_D = 0.01;
rho = 1.225; % [kg/m^3] air density
g = 9.81; % [m/s^2]
As = pi * (0.07874/2)^2; % [m^2] cross-sectional area
tb = 1.2; % [s] burnout time, pulled from given plot of thrust vs. time:
% https://wildmanrocketry.com/products/g74-6w?_pos=3&_sid=99a5eab6f&_ss=r
t_chute = tb + 6; % [s] time chute deploys after burnout

t = [0]; dt = 0.01;
v = [0];
theta = [pi/6];
x = [0];
h = [0];
L = [0];
imax = 100000;

for i = 1:imax
    if t(i) < tb
        F = 78;
        m = 0.921 - 0.087 * t(i);
    else
        F = 0;
    end

    Vrel = sqrt( (v(i)*sin(theta(i)))^2 + (Vw + v(i)*cos(theta(i)))^2 );
    D = 0.5 * rho * Vrel^2 * C_D * As;
    L(i) = 0.5 * rho * Vrel^2 * C_L * As;
    psi = atan( v(i)*sin(theta(i)) / (Vw + v(i)*cos(theta(i))) );
    G = ((F - D)*cos(psi - theta(i)) - L(i)*sin(psi - theta(i)) - ...
         m*g*sin(theta(i))) / m;
    dist = sqrt( x(i)^2 + h(i)^2 );
    if dist < Lrod
        H = 0;
    else
        H = ((F - D)*sin(psi - theta(i)) + L(i)*cos(psi - theta(i)) - ...
             m*g*cos(theta(i))) / (m * v(i));
    end

    t(i+1) = t(i) + dt;
    if t(i+1) < tb
        Fstar = 78;
        mstar = 0.921 - 0.087 * t(i+1);
    else
        Fstar = 0;
    end

    vstar = v(i) + G*dt;
    thetastar = theta(i) + H*dt;
    Vrelstar = sqrt( (vstar*sin(thetastar))^2 + (Vw + vstar*cos(thetastar))^2 );
    Dstar = 0.5 * rho * Vrelstar^2 * C_D * As;
    Lstar = 0.5 * rho * Vrelstar^2 * C_L * As;
    psistar = atan( vstar*sin(thetastar) / (Vw + vstar*cos(thetastar)) );
    Gstar = ((Fstar - Dstar)*cos(psistar - thetastar) - Lstar*sin(psistar - thetastar) - ...
         mstar*g*sin(thetastar)) / mstar;
    if dist < Lrod
        Hstar = 0;
    else
        Hstar = ((Fstar - Dstar)*sin(psistar - thetastar) + Lstar*cos(psistar - thetastar) - ...
             mstar*g*cos(thetastar)) / (mstar * vstar);
    end

    v(i+1) = v(i) + 0.5*(G+Gstar)*dt;
    theta(i+1) = theta(i) + 0.5*(H+Hstar)*dt;
    x(i+1) = x(i) + 0.5*( v(i)*cos(theta(i)) + v(i+1)*cos(theta(i+1)) )*dt;
    h(i+1) = h(i) + 0.5*( v(i)*sin(theta(i)) + v(i+1)*sin(theta(i+1)) )*dt;

    if h(i+1) < 0
        i = imax; %#ok<FXSET>
    elseif t(i) > 7.5
        break
    end

end

t = t.';
v = v.';
theta = theta.';
x = x.';
h = h.';