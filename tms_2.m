
clear; clc; close all;

%% License and Environment Check
matlab_ver = ver;
simscape_licensed = license('test', 'Simscape');
simulink_licensed = license('test', 'Simulink');
fprintf('MATLAB Version: %s\n', matlab_ver(1).Version);
fprintf('Simscape Licensed: %s\n', mat2str(simscape_licensed));
fprintf('Simulink Licensed: %s\n', mat2str(simulink_licensed));

% Check Simulink and Integrator block
use_simulink = simulink_licensed;
if use_simulink
    try
        add_block('simulink/Continuous/Integrator', 'test/Integrator');
        delete_block('test/Integrator');
        bdclose('test');
        fprintf('Simulink Integrator block available.\n');
    catch e
        warning('Simulink Integrator block unavailable: %s. Using numerical solver.', e.message);
        use_simulink = false;
    end
else
    fprintf('Simulink not licensed. Using numerical ODE solver.\n');
end

%% Battery Pack Parameters
% Cell specs (21700 lithium-ion)
cell_voltage = 3.6; % Nominal voltage (V)
cell_capacity = 8.55; % Capacity (Ah, adjusted for ~4kWh)
cell_resistance_base = 0.007; % Base internal resistance (Ohm, 7 mOhm)
cell_mass = 0.07; % Mass per cell (kg)
cell_specific_heat = 900; % Specific heat (J/kg-K)

% Pack configuration: 13s10p
num_series = 13;
num_parallel = 10;
pack_voltage = num_series * cell_voltage; % 46.8V
pack_capacity = num_parallel * cell_capacity; % 85.5Ah
pack_energy = pack_voltage * pack_capacity / 1000; % ~4kWh

% Cell resistance variation (±10%)
rng(42); % For reproducibility
cell_resistance = cell_resistance_base * (1 + 0.1 * (2 * rand(num_series, num_parallel) - 1));
avg_cell_resistance = mean(cell_resistance(:)); % Average for heat calculation

% Ambient conditions
T_ambient_base = 35; % °C (nominal ambient)
T_ambient_worst = 40; % °C (worst-case test)
T_max_limit = 45; % °C
T_runaway = 60; % °C (thermal runaway threshold)
T_coolant_in = 24; % °C (water-glycol, lowered for stronger cooling)
T_initial = 25; % °C (room temperature)

%% Electrical Load Profile (Drive Cycle Approximation)
% Piecewise constant C-rate: 0.5C, 2C, 3C, shorter segments
time_segments = [0 150 300 450 600 750 900 1050 1200 1350 1500 1650 1800 1950 2100 2250 2400 2550 2700 2850 3000 3150 3300 3450 3600];
c_rates = [0.5 2 1 3 0.5 2 1 3 0.5 2 1 3 0.5 2 1 3 0.5 2 1 3 0.5 2 1 3 0.5]; % Aggressive cycle
current_total = c_rates * pack_capacity; % Total current (A)
current_cell = current_total / num_parallel; % Per cell (A)

% Precompute heat generation per segment
heat_gen_cell = zeros(1, length(c_rates));
for i = 1:length(c_rates)
    heat_gen_cell(i) = (current_cell(i)^2) * avg_cell_resistance + 0.5; % I²R + entropic
end
total_heat = heat_gen_cell * num_series * num_parallel; % Total heat (W)

%% Cooling System Parameters
% Cooling plate (50/50 water-glycol mix)
pipe_surface_area = 0.2; % m² (multi-pass plate)
h_conv_base = 1000; % Base convective coefficient (W/m²-K, increased for better cooling)
R_thermal = 0.02; % Cell-to-cold-plate thermal resistance (K/W, reduced)
tau_delay = 10; % Coolant delay time constant (s, reduced for faster response)
m_dot = 0.02; % Assumed coolant mass flow (kg/s, increased for DeltaT)
Cp_coolant = 3500; % Specific heat of water-glycol (J/kg-K)

% Practical additions
mass_factor = 1.2; % 20% increase for casing and coolant mass
ambient_noise_std = 1.0; % Increased ambient temperature noise (±1°C)

% PID control (tuned for tighter ambient tracking)
k_P = 120; % Proportional gain (W/m²-K/°C, increased)
k_I = 1.5; % Integral gain (W/m²-K/°C/s, increased for steady-state error)
k_D = 25; % Derivative gain (W/m²-K·s/°C, increased)

%% Simulation (Nominal Case: T_ambient_base = 35°C)
mCp = mass_factor * cell_mass * num_series * num_parallel * cell_specific_heat; % J/K
T_max = NaN; T_final = NaN; T_coolant_out = NaN;
temp_profile = []; time_data = []; runaway_flag = false;

if use_simulink
    % Create Simulink model
    model_name = 'BatteryCoolingModel';
    bdclose(model_name);
    new_system(model_name);
    % open_system(model_name); % Uncomment for debugging

    % Blocks
    add_block('simulink/Continuous/Integrator', [model_name '/ThermalDynamics']);
    set_param([model_name '/ThermalDynamics'], 'InitialCondition', num2str(T_initial));
    add_block('simulink/Sources/Signal Builder', [model_name '/HeatSource']);
    set_param([model_name '/HeatSource'], 'Time', sprintf('[%s]', num2str(time_segments)), ...
        'Signals', sprintf('[%s]', num2str(total_heat)));
    add_block('simulink/Sources/Random Number', [model_name '/AmbientNoise']);
    set_param([model_name '/AmbientNoise'], 'Mean', '0', 'Variance', num2str(ambient_noise_std^2), 'SampleTime', '1');
    add_block('simulink/Sources/Constant', [model_name '/AmbientBase']);
    set_param([model_name '/AmbientBase'], 'Value', num2str(T_ambient_base));
    add_block('simulink/Commonly Used Blocks/Sum', [model_name '/SumAmbient']);
    set_param([model_name '/SumAmbient'], 'Inputs', '++');
    add_block('simulink/Continuous/PID Controller', [model_name '/PIDController']);
    set_param([model_name '/PIDController'], ...
        'P', num2str(k_P * pipe_surface_area), ...
        'I', num2str(k_I * pipe_surface_area), ...
        'D', num2str(k_D * pipe_surface_area));
    add_block('simulink/Commonly Used Blocks/Sum', [model_name '/TempError']);
    set_param([model_name '/TempError'], 'Inputs', '+-');
    add_block('simulink/Continuous/Transfer Fcn', [model_name '/CoolantDelay']);
    set_param([model_name '/CoolantDelay'], 'Numerator', '[1]', 'Denominator', sprintf('[%d 1]', tau_delay));
    add_block('simulink/Commonly Used Blocks/Gain', [model_name '/BaseCoolingGain']);
    set_param([model_name '/BaseCoolingGain'], 'Gain', num2str(-pipe_surface_area));
    add_block('simulink/Sources/Constant', [model_name '/CoolantTemp']);
    set_param([model_name '/CoolantTemp'], 'Value', num2str(T_coolant_in));
    add_block('simulink/Commonly Used Blocks/Sum', [model_name '/SumBaseCooling']);
    set_param([model_name '/SumBaseCooling'], 'Inputs', '+-');
    add_block('simulink/Math Operations/Math Function', [model_name '/NonlinearCooling']);
    set_param([model_name '/NonlinearCooling'], 'Function', 'recip');
    add_block('simulink/Math Operations/Sum', [model_name '/NonlinearFactor']);
    set_param([model_name '/NonlinearFactor'], 'Inputs', '++');
    add_block('simulink/Sources/Constant', [model_name '/NonlinearBase']);
    set_param([model_name '/NonlinearBase'], 'Value', '1');
    add_block('simulink/Commonly Used Blocks/Gain', [model_name '/NonlinearGain']);
    set_param([model_name '/NonlinearGain'], 'Gain', '0.01');
    add_block('simulink/Commonly Used Blocks/Product', [model_name '/ApplyNonlinear']);
    set_param([model_name '/ApplyNonlinear'], 'Inputs', '2');
    add_block('simulink/Commonly Used Blocks/Sum', [model_name '/TotalCooling']);
    set_param([model_name '/TotalCooling'], 'Inputs', '+-');
    add_block('simulink/Commonly Used Blocks/Gain', [model_name '/ThermalResistance']);
    set_param([model_name '/ThermalResistance'], 'Gain', num2str(1 / (1/h_conv_base + R_thermal)));
    add_block('simulink/Commonly Used Blocks/Gain', [model_name '/ThermalMass']);
    set_param([model_name '/ThermalMass'], 'Gain', num2str(1/mCp));
    add_block('simulink/Sinks/Scope', [model_name '/TempScope']);
    set_param([model_name '/TempScope'], 'SaveToWorkspace', 'on', 'VariableName', 'TempData');

    % Connections
    add_line(model_name, 'AmbientBase/1', 'SumAmbient/1', 'autorouting', 'on');
    add_line(model_name, 'AmbientNoise/1', 'SumAmbient/2', 'autorouting', 'on');
    add_line(model_name, 'SumAmbient/1', 'TempError/2', 'autorouting', 'on');
    add_line(model_name, 'ThermalDynamics/1', 'TempError/1', 'autorouting', 'on');
    add_line(model_name, 'TempError/1', 'PIDController/1', 'autorouting', 'on');
    add_line(model_name, 'PIDController/1', 'CoolantDelay/1', 'autorouting', 'on');
    add_line(model_name, 'ThermalDynamics/1', 'BaseCoolingGain/1', 'autorouting', 'on');
    add_line(model_name, 'BaseCoolingGain/1', 'SumBaseCooling/1', 'autorouting', 'on');
    add_line(model_name, 'CoolantTemp/1', 'SumBaseCooling/2', 'autorouting', 'on');
    add_line(model_name, 'ThermalDynamics/1', 'NonlinearGain/1', 'autorouting', 'on');
    add_line(model_name, 'NonlinearGain/1', 'NonlinearFactor/1', 'autorouting', 'on');
    add_line(model_name, 'NonlinearBase/1', 'NonlinearFactor/2', 'autorouting', 'on');
    add_line(model_name, 'NonlinearFactor/1', 'NonlinearCooling/1', 'autorouting', 'on');
    add_line(model_name, 'SumBaseCooling/1', 'ApplyNonlinear/1', 'autorouting', 'on');
    add_line(model_name, 'NonlinearCooling/1', 'ApplyNonlinear/2', 'autorouting', 'on');
    add_line(model_name, 'ApplyNonlinear/1', 'ThermalResistance/1', 'autorouting', 'on');
    add_line(model_name, 'ThermalResistance/1', 'TotalCooling/2', 'autorouting', 'on');
    add_line(model_name, 'CoolantDelay/1', 'TotalCooling/1', 'autorouting', 'on');
    add_line(model_name, 'HeatSource/1', 'ThermalMass/1', 'autorouting', 'on');
    add_line(model_name, 'TotalCooling/1', 'ThermalMass/2', 'autorouting', 'on');
    add_line(model_name, 'ThermalMass/1', 'ThermalDynamics/1', 'autorouting', 'on');
    add_line(model_name, 'ThermalDynamics/1', 'TempScope/1', 'autorouting', 'on');

    % Simulation settings
    set_param(model_name, ...
        'StopTime', '3600', ...
        'Solver', 'ode23t', ...
        'RelTol', '1e-4', ...
        'AbsTol', '1e-6');

    % Simulate (Nominal Case)
    try
        sim(model_name);
        temp_var = evalin('base', 'TempData');
        T_max = max(temp_var.signals.values);
        T_final = temp_var.signals.values(end);
        temp_profile = temp_var.signals.values;
        time_data = temp_var.time;
        runaway_flag = any(temp_var.signals.values > T_runaway);
        % Estimate coolant DeltaT
        Q_cool_avg = mean(h_conv_base * pipe_surface_area * (temp_var.signals.values - T_coolant_in));
        T_coolant_out = T_coolant_in + Q_cool_avg / (m_dot * Cp_coolant);
    catch e
        fprintf('Nominal simulation failed: %s\n', e.message);
        T_max = NaN; T_final = NaN; temp_profile = []; time_data = [];
        bdclose(model_name);
    end

    bdclose(model_name);
else
    % Numerical ODE solver (RK4)
    t = 0:0.01:3600; dt = 0.01;
    T = zeros(size(t)); T(1) = T_initial;
    Q_cool_delayed = 0; integral_error = 0; prev_error = 0;
    T_coolant_out_sum = 0;

    for j = 2:length(t)
        % Interpolate heat generation
        total_heat_t = interp1(time_segments, total_heat, t(j-1), 'previous');
        % Ambient temperature with noise
        T_ambient = T_ambient_base + ambient_noise_std * randn();
        % RK4 steps
        k1 = thermal_dynamics(T(j-1), t(j-1), total_heat_t, T_ambient, T_coolant_in, ...
            h_conv_base, pipe_surface_area, R_thermal, k_P, k_I, k_D, ...
            integral_error, prev_error, Q_cool_delayed, dt, tau_delay, mCp);
        k2 = thermal_dynamics(T(j-1) + dt/2*k1, t(j-1) + dt/2, total_heat_t, T_ambient, ...
            T_coolant_in, h_conv_base, pipe_surface_area, R_thermal, k_P, k_I, k_D, ...
            integral_error, prev_error, Q_cool_delayed, dt, tau_delay, mCp);
        k3 = thermal_dynamics(T(j-1) + dt/2*k2, t(j-1) + dt/2, total_heat_t, T_ambient, ...
            T_coolant_in, h_conv_base, pipe_surface_area, R_thermal, k_P, k_I, k_D, ...
            integral_error, prev_error, Q_cool_delayed, dt, tau_delay, mCp);
        k4 = thermal_dynamics(T(j-1) + dt*k3, t(j-1) + dt, total_heat_t, T_ambient, ...
            T_coolant_in, h_conv_base, pipe_surface_area, R_thermal, k_P, k_I, k_D, ...
            integral_error, prev_error, Q_cool_delayed, dt, tau_delay, mCp);
        % Update temperature
        T(j) = T(j-1) + (dt/6) * (k1 + 2*k2 + 2*k3 + k4);
        % Update PID and cooling
        error = T(j-1) - T_ambient;
        integral_error = integral_error + error * dt;
        derivative_error = (error - prev_error) / dt;
        h_control = k_P * error + k_I * integral_error + k_D * derivative_error;
        prev_error = error;
        delta_T = T(j-1) - T_coolant_in;
        h_eff = h_conv_base / (1 + 0.01 * delta_T);
        h_total = h_eff + max(0, h_control);
        Q_cool = (T(j-1) - T_coolant_in) / (R_thermal + 1/(h_total * pipe_surface_area));
        Q_cool_delayed = Q_cool_delayed + (Q_cool - Q_cool_delayed) * (dt / tau_delay);
        % Coolant DeltaT
        T_coolant_out_sum = T_coolant_out_sum + Q_cool;
        % Runaway check
        if T(j) > T_runaway
            runaway_flag = true;
        end
    end
    T_max = max(T);
    T_final = T(end);
    temp_profile = T;
    time_data = t;
    T_coolant_out = T_coolant_in + (T_coolant_out_sum / length(t)) / (m_dot * Cp_coolant);
end

%% Simulation (Worst-Case: T_ambient_base = 40°C)
T_max_worst = NaN; T_final_worst = NaN; runaway_flag_worst = false;
if use_simulink
    % Reuse model, update ambient temperature
    model_name = 'BatteryCoolingModel';
    bdclose(model_name);
    new_system(model_name);
    add_block('simulink/Continuous/Integrator', [model_name '/ThermalDynamics']);
    set_param([model_name '/ThermalDynamics'], 'InitialCondition', num2str(T_initial));
    add_block('simulink/Sources/Signal Builder', [model_name '/HeatSource']);
    set_param([model_name '/HeatSource'], 'Time', sprintf('[%s]', num2str(time_segments)), ...
        'Signals', sprintf('[%s]', num2str(total_heat)));
    add_block('simulink/Sources/Random Number', [model_name '/AmbientNoise']);
    set_param([model_name '/AmbientNoise'], 'Mean', '0', 'Variance', num2str(ambient_noise_std^2), 'SampleTime', '1');
    add_block('simulink/Sources/Constant', [model_name '/AmbientBase']);
    set_param([model_name '/AmbientBase'], 'Value', num2str(T_ambient_worst));
    add_block('simulink/Commonly Used Blocks/Sum', [model_name '/SumAmbient']);
    set_param([model_name '/SumAmbient'], 'Inputs', '++');
    add_block('simulink/Continuous/PID Controller', [model_name '/PIDController']);
    set_param([model_name '/PIDController'], ...
        'P', num2str(k_P * pipe_surface_area), ...
        'I', num2str(k_I * pipe_surface_area), ...
        'D', num2str(k_D * pipe_surface_area));
    add_block('simulink/Commonly Used Blocks/Sum', [model_name '/TempError']);
    set_param([model_name '/TempError'], 'Inputs', '+-');
    add_block('simulink/Continuous/Transfer Fcn', [model_name '/CoolantDelay']);
    set_param([model_name '/CoolantDelay'], 'Numerator', '[1]', 'Denominator', sprintf('[%d 1]', tau_delay));
    add_block('simulink/Commonly Used Blocks/Gain', [model_name '/BaseCoolingGain']);
    set_param([model_name '/BaseCoolingGain'], 'Gain', num2str(-pipe_surface_area));
    add_block('simulink/Sources/Constant', [model_name '/CoolantTemp']);
    set_param([model_name '/CoolantTemp'], 'Value', num2str(T_coolant_in));
    add_block('simulink/Commonly Used Blocks/Sum', [model_name '/SumBaseCooling']);
    set_param([model_name '/SumBaseCooling'], 'Inputs', '+-');
    add_block('simulink/Math Operations/Math Function', [model_name '/NonlinearCooling']);
    set_param([model_name '/NonlinearCooling'], 'Function', 'recip');
    add_block('simulink/Math Operations/Sum', [model_name '/NonlinearFactor']);
    set_param([model_name '/NonlinearFactor'], 'Inputs', '++');
    add_block('simulink/Sources/Constant', [model_name '/NonlinearBase']);
    set_param([model_name '/NonlinearBase'], 'Value', '1');
    add_block('simulink/Commonly Used Blocks/Gain', [model_name '/NonlinearGain']);
    set_param([model_name '/NonlinearGain'], 'Gain', '0.01');
    add_block('simulink/Commonly Used Blocks/Product', [model_name '/ApplyNonlinear']);
    set_param([model_name '/ApplyNonlinear'], 'Inputs', '2');
    add_block('simulink/Commonly Used Blocks/Sum', [model_name '/TotalCooling']);
    set_param([model_name '/TotalCooling'], 'Inputs', '+-');
    add_block('simulink/Commonly Used Blocks/Gain', [model_name '/ThermalResistance']);
    set_param([model_name '/ThermalResistance'], 'Gain', num2str(1 / (1/h_conv_base + R_thermal)));
    add_block('simulink/Commonly Used Blocks/Gain', [model_name '/ThermalMass']);
    set_param([model_name '/ThermalMass'], 'Gain', num2str(1/mCp));
    add_block('simulink/Sinks/Scope', [model_name '/TempScope']);
    set_param([model_name '/TempScope'], 'SaveToWorkspace', 'on', 'VariableName', 'TempDataWorst');

    % Connections
    add_line(model_name, 'AmbientBase/1', 'SumAmbient/1', 'autorouting', 'on');
    add_line(model_name, 'AmbientNoise/1', 'SumAmbient/2', 'autorouting', 'on');
    add_line(model_name, 'SumAmbient/1', 'TempError/2', 'autorouting', 'on');
    add_line(model_name, 'ThermalDynamics/1', 'TempError/1', 'autorouting', 'on');
    add_line(model_name, 'TempError/1', 'PIDController/1', 'autorouting', 'on');
    add_line(model_name, 'PIDController/1', 'CoolantDelay/1', 'autorouting', 'on');
    add_line(model_name, 'ThermalDynamics/1', 'BaseCoolingGain/1', 'autorouting', 'on');
    add_line(model_name, 'BaseCoolingGain/1', 'SumBaseCooling/1', 'autorouting', 'on');
    add_line(model_name, 'CoolantTemp/1', 'SumBaseCooling/2', 'autorouting', 'on');
    add_line(model_name, 'ThermalDynamics/1', 'NonlinearGain/1', 'autorouting', 'on');
    add_line(model_name, 'NonlinearGain/1', 'NonlinearFactor/1', 'autorouting', 'on');
    add_line(model_name, 'NonlinearBase/1', 'NonlinearFactor/2', 'autorouting', 'on');
    add_line(model_name, 'NonlinearFactor/1', 'NonlinearCooling/1', 'autorouting', 'on');
    add_line(model_name, 'SumBaseCooling/1', 'ApplyNonlinear/1', 'autorouting', 'on');
    add_line(model_name, 'NonlinearCooling/1', 'ApplyNonlinear/2', 'autorouting', 'on');
    add_line(model_name, 'ApplyNonlinear/1', 'ThermalResistance/1', 'autorouting', 'on');
    add_line(model_name, 'ThermalResistance/1', 'TotalCooling/2', 'autorouting', 'on');
    add_line(model_name, 'CoolantDelay/1', 'TotalCooling/1', 'autorouting', 'on');
    add_line(model_name, 'HeatSource/1', 'ThermalMass/1', 'autorouting', 'on');
    add_line(model_name, 'TotalCooling/1', 'ThermalMass/2', 'autorouting', 'on');
    add_line(model_name, 'ThermalMass/1', 'ThermalDynamics/1', 'autorouting', 'on');
    add_line(model_name, 'ThermalDynamics/1', 'TempScope/1', 'autorouting', 'on');

    % Simulate
    try
        sim(model_name);
        temp_var_worst = evalin('base', 'TempDataWorst');
        T_max_worst = max(temp_var_worst.signals.values);
        T_final_worst = temp_var_worst.signals.values(end);
        runaway_flag_worst = any(temp_var_worst.signals.values > T_runaway);
    catch e
        fprintf('Worst-case simulation failed: %s\n', e.message);
        T_max_worst = NaN; T_final_worst = NaN;
        bdclose(model_name);
    end

    bdclose(model_name);
else
    % Numerical ODE solver (Worst-Case, RK4)
    t = 0:0.01:3600; dt = 0.01;
    T = zeros(size(t)); T(1) = T_initial;
    Q_cool_delayed = 0; integral_error = 0; prev_error = 0;
    for j = 2:length(t)
        total_heat_t = interp1(time_segments, total_heat, t(j-1), 'previous');
        T_ambient = T_ambient_worst + ambient_noise_std * randn();
        k1 = thermal_dynamics(T(j-1), t(j-1), total_heat_t, T_ambient, T_coolant_in, ...
            h_conv_base, pipe_surface_area, R_thermal, k_P, k_I, k_D, ...
            integral_error, prev_error, Q_cool_delayed, dt, tau_delay, mCp);
        k2 = thermal_dynamics(T(j-1) + dt/2*k1, t(j-1) + dt/2, total_heat_t, T_ambient, ...
            T_coolant_in, h_conv_base, pipe_surface_area, R_thermal, k_P, k_I, k_D, ...
            integral_error, prev_error, Q_cool_delayed, dt, tau_delay, mCp);
        k3 = thermal_dynamics(T(j-1) + dt/2*k2, t(j-1) + dt/2, total_heat_t, T_ambient, ...
            T_coolant_in, h_conv_base, pipe_surface_area, R_thermal, k_P, k_I, k_D, ...
            integral_error, prev_error, Q_cool_delayed, dt, tau_delay, mCp);
        k4 = thermal_dynamics(T(j-1) + dt*k3, t(j-1) + dt, total_heat_t, T_ambient, ...
            T_coolant_in, h_conv_base, pipe_surface_area, R_thermal, k_P, k_I, k_D, ...
            integral_error, prev_error, Q_cool_delayed, dt, tau_delay, mCp);
        T(j) = T(j-1) + (dt/6) * (k1 + 2*k2 + 2*k3 + k4);
        error = T(j-1) - T_ambient;
        integral_error = integral_error + error * dt;
        derivative_error = (error - prev_error) / dt;
        h_control = k_P * error + k_I * integral_error + k_D * derivative_error;
        prev_error = error;
        delta_T = T(j-1) - T_coolant_in;
        h_eff = h_conv_base / (1 + 0.01 * delta_T);
        h_total = h_eff + max(0, h_control);
        Q_cool = (T(j-1) - T_coolant_in) / (R_thermal + 1/(h_total * pipe_surface_area));
        Q_cool_delayed = Q_cool_delayed + (Q_cool - Q_cool_delayed) * (dt / tau_delay);
        if T(j) > T_runaway
            runaway_flag_worst = true;
        end
    end
    T_max_worst = max(T);
    T_final_worst = T(end);
end

%% Results (Nominal Case)
fprintf('\nNominal Case (Ambient = %.0f°C):\n', T_ambient_base);
fprintf('Cell Temperature (Max): %.2f °C\n', T_max);
fprintf('Cell Temperature (Final): %.2f °C\n', T_final);
fprintf('Coolant Inlet Temperature: %.2f °C\n', T_coolant_in);
fprintf('Coolant Outlet Temperature (Avg): %.2f °C\n', T_coolant_out);
fprintf('Ambient Temperature (Nominal): %.2f °C\n', T_ambient_base);
if T_max <= T_max_limit && ~isnan(T_max)
    fprintf('Cell temperature is within limit (%.0f°C).\n', T_max_limit);
    if abs(T_final - T_ambient_base) <= 1
        fprintf('Cell temperature stabilized within ambient range (%.0f ± 1°C).\n', T_ambient_base);
    else
        fprintf('Cell temperature stabilized outside ambient range (%.0f ± 1°C) but within limit.\n', T_ambient_base);
    end
else
    warning('Cell temperature exceeds limit (%.0f°C) or simulation failed.', T_max_limit);
end
if runaway_flag
    warning('Thermal runaway risk detected: Temperature exceeded %.0f°C.', T_runaway);
end

%% Results (Worst-Case)
fprintf('\nWorst-Case (Ambient = %.0f°C):\n', T_ambient_worst);
fprintf('Cell Temperature (Max): %.2f °C\n', T_max_worst);
fprintf('Cell Temperature (Final): %.2f °C\n', T_final_worst);
fprintf('Coolant Inlet Temperature: %.2f °C\n', T_coolant_in);
fprintf('Ambient Temperature (Nominal): %.2f °C\n', T_ambient_worst);
if T_max_worst <= T_max_limit && ~isnan(T_max_worst)
    fprintf('Cell temperature is within limit (%.0f°C).\n', T_max_limit);
    if abs(T_final_worst - T_ambient_worst) <= 1
        fprintf('Cell temperature stabilized within ambient range (%.0f ± 1°C).\n', T_ambient_worst);
    else
        fprintf('Cell temperature stabilized outside ambient range (%.0f ± 1°C) but within limit.\n', T_ambient_worst);
    end
else
    warning('Cell temperature exceeds limit (%.0f°C) or simulation failed.', T_max_limit);
end
if runaway_flag_worst
    warning('Thermal runaway risk detected: Temperature exceeded %.0f°C.', T_runaway);
end

%% Plot: Transient Temperature Profile (Nominal Case)
figure('Name', 'Temperature Profile');
if ~isempty(temp_profile)
    plot(time_data, temp_profile, 'b-', 'LineWidth', 2, 'DisplayName', 'Cell Temp (Nominal)');
    hold on;
end
plot([0 3600], [T_coolant_in T_coolant_in], 'c--', 'LineWidth', 1.5, ...
    'DisplayName', sprintf('Coolant Inlet (%.0f°C)', T_coolant_in));
plot([0 3600], [T_ambient_base T_ambient_base], 'm--', 'LineWidth', 1.5, ...
    'DisplayName', sprintf('Ambient (%.0f°C)', T_ambient_base));
plot([0 3600], [45 45], 'r--', 'LineWidth', 1.5, 'DisplayName', '45°C Limit');
plot([0 3600], [60 60], 'k:', 'LineWidth', 1.5, 'DisplayName', 'Runaway Threshold (60°C)');
xlabel('Time (s)');
ylabel('Temperature (°C)');
title('Cell Temperature Over Time (Practical Model, Nominal Case)');
grid on;
legend('Location', 'best');

%% Save Results
save('BatteryCoolingResults.mat', 'T_max', 'T_final', 'T_coolant_in', 'T_coolant_out', ...
    'T_ambient_base', 'temp_profile', 'time_data', 'T_max_worst', 'T_final_worst', 'T_ambient_worst');

%% Local Function for RK4
function dT = thermal_dynamics(T, t, total_heat_t, T_ambient, T_coolant_in, ...
        h_conv_base, pipe_surface_area, R_thermal, k_P, k_I, k_D, ...
        integral_error, prev_error, Q_cool_delayed, dt, tau_delay, mCp)
    error = T - T_ambient;
    integral_error = integral_error + error * dt;
    derivative_error = (error - prev_error) / dt;
    h_control = k_P * error + k_I * integral_error + k_D * derivative_error;
    delta_T = T - T_coolant_in;
    h_eff = h_conv_base / (1 + 0.01 * delta_T);
    h_total = h_eff + max(0, h_control);
    Q_cool = (T - T_coolant_in) / (R_thermal + 1/(h_total * pipe_surface_area));
    Q_cool_delayed = Q_cool_delayed + (Q_cool - Q_cool_delayed) * (dt / tau_delay);
    dT = (total_heat_t - Q_cool_delayed) / mCp;
end

