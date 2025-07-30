function [t,y1,y2,y3,y4,y5] = KM_solver(PR,equilb,sig,mu_l,c_l,gam_v,p_init,p_v,R0,grow,rho_l,t0)
% KM_solver - Solves the Keller-Miksis equation for bubble dynamics
%
% Inputs:
%   PR      : Pressure ratio (p_inf / 1e5)
%   equilb  : Equilibrium flag (0: R_dot = 0, 1: initial R_dot estimated)
%   sig     : Surface tension [N/m]
%   mu_l    : Liquid viscosity [Pa·s]
%   c_l     : Speed of sound in the liquid [m/s]
%   gam_v   : Specific heat ratio of gas inside the bubble
%   p_init  : Initial bubble pressure [Pa]
%   p_v     : Vapor pressure [Pa]
%   R0      : Initial bubble radius [m]
%   grow    : Bubble growth or collapse flag (0: growth, 1: collapse)
%   rho_l   : Liquid density [kg/m^3]
%   t0      : Initial time [s]
%
% Outputs:
%   t   : Time vector [s]
%   y1  : Bubble radius [m]
%   y2  : Wall velocity [m/s]
%   y3  : Bubble pressure [Pa]
%   y4  : Bubble wall acceleration [m/s²]
%   y5  : Time derivative of bubble pressure [Pa/s]

%% Derived parameters
p_inf = PR * 1e5;  % Far-field pressure [Pa]

% Determine initial bubble wall velocity based on equilibrium assumption
if equilb == 0
    R_dot = 0;
else
    R_dot = (p_init - p_inf) / (rho_l * c_l);
end

%% Initial bubble pressure
pg = p_init;     % Gas pressure
pv = p_v;        % Vapor pressure
pb = pg + pv;    % Total initial bubble pressure

%% Characteristic scales for non-dimensionalization
Lc = R0;                            % Characteristic length
tc = Lc * sqrt(rho_l / p_inf);      % Characteristic time
Vc = Lc / tc;                       % Characteristic velocity
p_inf_st = p_inf / (rho_l * Vc^2);  % Non-dimensional far-field pressure

% Convert initial time and compute final simulation time
t_init = t0 / tc;
t_final = t0 + 3.5 * tc;

%% Dimensionless parameters
Re_l = (rho_l * Vc * Lc) / mu_l;   % Reynolds number
We   = (rho_l * Vc^2 * Lc) / sig;  % Weber number
c_st = c_l / Vc;                   % Non-dimensional sound speed

%% Output arrays for tracking variables
t_tmp = 0.0; 
tout = []; 
pinfout = []; 
Twout = []; 
k1_wout = [];
R_2dot_save = []; 
p_b_dot_save = [];

%% Initial conditions (non-dimensional)
y0(1,1) = R0 / R0;              % R/R0 (starts as 1)
y0(1,2) = R_dot / Vc;           % Non-dimensional wall velocity
y0(1,3) = pb / (rho_l * Vc^2);  % Non-dimensional bubble pressure

%% Time integration using ode15s
tspan = [t_init t_final / tc];
options = odeset('RelTol',1e-16,'AbsTol',1e-16);
[t, y] = ode15s(@derivs, tspan, y0, options);  % Solve ODE system

%% Rescale outputs to physical units
t  = t * tc;
y1 = y(:,1) * Lc;                 % Radius [m]
y2 = y(:,2) * Vc;                 % Wall velocity [m/s]
y3 = y(:,3) * rho_l * Vc^2;       % Bubble pressure [Pa]
y4 = R_2dot_save;                 % Wall acceleration [m/s²]
y5 = p_b_dot_save;                % Pressure derivative [Pa/s]

%% Nested ODE function
function f = derivs(t, y)
    R_st   = y(1);      % Non-dimensional radius
    R_dot  = y(2);      % Non-dimensional wall velocity
    p_b_st = y(3);      % Non-dimensional pressure inside bubble

    % Compute time derivative of bubble pressure
    f(3,1) = -3.0 / R_st * (gam_v * p_b_st * R_dot);
    p_b_dot = f(3,1);

    % R_dot is just the velocity
    f(1,1) = R_dot;

    % Keller-Miksis equation for R_2dot
    f(2,1) = (-3/2 * (1 - R_dot / (3 * c_st)) * R_dot^2 ...
             + R_st / c_st * p_b_dot ...
             - 4 / Re_l * R_dot / R_st ...
             - 2 / We / R_st ...
             + (p_b_st - p_inf_st) * (1 + R_dot / c_st)) ...
             / (R_st * (1 - R_dot / c_st) + 4 / Re_l / c_st);

    R_2dot = f(2,1);  % Save acceleration for output

    % Avoid duplicate time entries in output
    Nt_tmp = numel(tout);
    if (Nt_tmp > 1) && (tout(Nt_tmp) == tout(Nt_tmp - 1))
        Nt_tmp = Nt_tmp - 1;
    end

    % Save output data at each time step
    if t ~= t_tmp
        if t <= t_tmp
            tout(Nt_tmp,1) = t;
            pinfout(Nt_tmp,1) = p_inf;
            R_2dot_save(Nt_tmp,1) = R_2dot;
            p_b_dot_save(Nt_tmp,1) = p_b_dot;
        else
            tout(Nt_tmp + 1,1) = t;
            pinfout(Nt_tmp + 1,1) = p_inf;
            R_2dot_save(Nt_tmp + 1,1) = R_2dot;
            p_b_dot_save(Nt_tmp + 1,1) = p_b_dot;
        end
    end

    t_tmp = t;  % Update previous time value
end

%% Utility functions (unused but defined for possible extensions)

% Estimate initial velocity and temperature based on pressure match
function [R_dot, T_b_init] = f_init()
    it = 0;
    while true
        it = it + 1;
        p_v_sat = f_pvsat(T_inf + dT_tmp);
        p_b = p_v_sat;
        R_dot_iner = sqrt((2/3) * (p_b - p_inf) / rho_l);
        if abs(R_dot_goal - R_dot_iner) < 0.1
            dT_super = dT_tmp;
            break;
        else
            dT_tmp = dT_tmp + 0.1;
        end
    end
    R_dot = R_dot_iner;
    T_b_init = T_inf + dT_tmp;
end

% Saturated vapor pressure [Pa] as a function of temperature [K]
function pvsattmp = f_pvsat(Ttmp)
    pvsattmp = 1.17e11 * exp(-5200.0 / Ttmp);
end

% Thermal conductivity of gas mixture (linear model)
function K1tmp = f_K1(Ttmp)
    K1tmp = Ak * Ttmp + Bk;
end

% Convert thermal time constant (tau) to temperature
function Ttmp = f_tautoT(tautmp)
    Ttmp = (sqrt(K1_0^2 + 2.0 * Ak * tautmp) - Bk) / Ak;
end

% Convert temperature to thermal time constant (tau)
function tautmp = f_Ttotau(Ttmp)
    tautmp = 0.5 * Ak * (Ttmp^2 - T0^2) + Bk * (Ttmp - T0);
end

end
