% run RPsolver
clear all; clc; close all; format compact;

%=====================================================================
% Run RPsolver for analyzing bubble dynamics
%=====================================================================
% Output variables:
% t  : time
% y1 : bubble radius
% y2 : bubble wall velocity
% y3 : bubble pressure

% Define pressure ratio input range
AA = [-1.2:0.2:2];
PRi = 10.^AA;

% Initial bubble radius [m]
R0  = 0.5/10^3;

% Select model settings
visurf = 0;     % Surface and viscosity model flag
                % 0: no viscosity, no surface tension
                % 1: with viscosity, no surface tension
                % 2: no viscosity, with surface tension
                % 3: with viscosity, with surface tension
equilb = 1;     % 0: equilibrium condition (R_dot = 0)
grow   = 1;     % 0: bubble growth, 1: bubble collapse

% Initial conditions and physical constants
p_init = 3550;     % Initial bubble pressure [Pa]
T_init = 300;      % Initial bubble temperature [K]
p_v    = 0;        % Vapor pressure [Pa]
t0     = 0;        % Initial time [s]
gam_v  = 1.4;      % Specific heat ratio of gas inside bubble
c_l    = 1510.0;   % Speed of sound in liquid [m/s]
rho_l  = 998;      % Liquid density [kg/m^3]

% Assign viscosity and surface tension based on model selection
if visurf == 0
    sig  = 1e-15;
    mu_l = 1e-15;
    name3 = 'no_pv_mu_sig\';
elseif visurf == 1
    sig  = 1e-15;
    mu_l = 8.3283e-4;
    name3 = 'no_pv_w_mu_n_sig\';
elseif visurf == 2
    sig  = 0.07197;
    mu_l = 1e-15;
    name3 = 'no_pv_n_mu_w_sig\';
elseif visurf == 3
    sig  = 0.07197;
    mu_l = 8.3283e-4;
    name3 = 'no_pv_w_mu_w_sig\';
end

% Loop over all pressure ratios
for i = 1:numel(PRi)
    PR     = PRi(i);    
    p_inf  = PR * 1e5;  % Far-field pressure [Pa]

    % Solve the Keller-Miksis (KM) equation
    [t,y1,y2,y3,y4,y5] = KM_solver(PR,equilb,sig,mu_l,c_l,gam_v,p_init,p_v,R0,grow,rho_l,t0);

    % Maximum bubble pressure and corresponding temperature
    pb_max(i,1) = max(y3);  
    Tb_max(i,1) = T_init * (p_init / pb_max(i))^((1 - gam_v) / gam_v);

    % Minimum bubble radius (as a multiple of R0) and corresponding time
    [Rmin, in] = min(y1/R0)
    Rmin_all(i,1) = Rmin;
    t(in);  % Display time at minimum radius

    % Plot setup (not executed here, only values assigned)
    t_plot = t * 1e6;     % Time in microseconds
    R_plot = y1 * 1e3;    % Radius in millimeters

    % Sample data to be saved (every third point)
    t_leng = numel(t_plot);
    out = [t_plot(1:3:t_leng) R_plot(1:3:t_leng)];

    % Define path for saving results based on flags
    if grow == 0 && equilb == 0
        name1 = ['C:\2_CAV\Local\0 R-P-type2\KM_c=const\1 grow_gas_equilb\', name3];
    elseif grow == 0 && equilb == 1
        name1 = ['C:\2_CAV\Local\0 R-P-type2\KM_c=const\1 grow_gas_non_equilb\', name3];
    elseif grow == 1 && equilb == 0
        name1 = ['C:\2_CAV\Local\0 R-P-type2\KM_c=const\2 col_gas_equilb\', name3];
    elseif grow == 1 && equilb == 1
        name1 = ['C:\2_CAV\Local\0 R-P-type2\KM_c=const\2 col_gas_non_equilb\', name3];
    end

    name2 = ['R0=', num2str(R0), '_PR=', num2str(PR), '\pb=', num2str(p_init/10^3), ...
             '_rho=', num2str(rho_l), '_gam=', num2str(gam_v)];
    mkdir([name1, name2]);  % Create output directory

    % % Save physical output variables
    % fileID0 = fopen([name1,name2,'\R_Q.txt'],'w');
    % fprintf(fileID0,'t(s)\t  R(m)\t  R_dot(m/s)\t  R_2dot(m/s^2)\t  p_b(pa)\t  p_b_dot(pa/s) \r\n');
    % 
    % fileID1 = fopen([name1,name2,'\R_Q_non.txt'],'w');
    % fprintf(fileID1,'t/t_non\t  t/t_c\t  R/R0\t  R_dot_non\t  R_2dot_non\t  p_b_non\t  p_b_dot_non\r\n');
    % 
    % output0 = [t y1 y2 y4 y3 y5];
    % fprintf(fileID0,'%.12e\t %.12e\t %.12e\t %.12e\t %.12e\t %.12e\r\n', output0');

    % Derive KM equation components
    R     = y1;
    R_dot = y2;
    p_l   = y3;

    % Estimate time derivative of pressure using gas dynamics
    p_l_dot = -3 * gam_v * R_dot ./ R .* p_l;

    % Compute second derivative of radius using KM equation (simplified & full form)
    R_2dot  = (-1.5*R_dot.^2 + (p_l - p_inf) / rho_l) ./ R;

    % Calculate forcing term for KM residual (diagnostic)
    F_2dot = R_dot.^3/2 + R .* R_dot .* R_2dot + R_dot .* (p_l - p_inf)/rho_l + R .* p_l_dot / rho_l;

    % KM equation components
    KM1 = R .* R_2dot;
    KM2 = 1.5 * R_dot.^2;
    KM3 = -F_2dot / c_l;
    KM4 = (p_l - p_inf) / rho_l;
    KM_tot = KM1 + KM2 + KM3 - KM4;

    % Diagnostic: check residuals of KM equation
    max(abs(KM_tot));
    min(abs(KM_tot));
    1;  % Marker (no-op)

    % Non-dimensionalization
    Lc   = R0;
    t_non = Lc * sqrt(rho_l / (p_inf - 0));  % Characteristic time
    tc    = 0.915 * t_non;
    Vc    = Lc / tc;

    R_non       = y1 / Lc;
    R_dot_non   = y2 / Vc;
    R_2dot_non  = y4 / (Lc / tc^2);
    p_b_non     = y3 / (rho_l * Vc^2);
    p_b_dot_non = y5 / (rho_l * Vc^2 / tc);

    output1 = [t/t_non t/tc R_non R_dot_non R_2dot_non p_b_non p_b_dot_non];
    % fprintf(fileID1,'%.8e\t %.8e\t %.8e\t %.8e\t %.8e\t %.8e\t %.8e\r\n', output0');
    % fclose(fileID0); 
    % fclose(fileID1);

    % Save max bubble pressure for each PR value
    pb_max_save(i) = max(y3);

end

% Compute corresponding maximum temperature for each PR
Tb_max_save = T_init * (p_init ./ pb_max_save).^((1 - gam_v) / gam_v);

% Prepare final output: [normalized PR, max pressure, max temperature]
pb_max_out = [(PRi' * 1e5 - 3550) / 3550, pb_max_save', Tb_max_save'];


