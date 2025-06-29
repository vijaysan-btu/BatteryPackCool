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
        add_block('simulink/Commonly Used Blocks/Integrator', 'test/Integrator');
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
cell_resistance = 0.02; % Internal resistance (Ohm)
cell_mass = 0.07; % Mass per cell (kg)
cell_specific_heat = 900; % Specific heat (J/kg-K)

% Pack configuration: 13s10p
num_series = 13;
num_parallel = 10;
pack_voltage = num_series * cell_voltage; % 46.8V
pack_capacity = num_parallel * cell_capacity; % 85.5Ah
pack_energy = pack_voltage * pack_capacity / 1000; % ~4kWh

% Discharge conditions
discharge_rate = 2; % 2C
current_total = discharge_rate * pack_capacity; % 171A
current_cell = current_total / num_parallel; % 17.1A per cell

% Heat generation per cell (I^2*R + entropic heating)
heat_gen_cell = (current_cell^2) * cell_resistance + 0.5; % ~5.84W + 0.5W entropic
total_heat_base = heat_gen_cell * num_series * num_parallel; % ~799.2W

% Ambient conditions
T_ambient_base = 35; % °C (nominal ambient)
T_max_limit = 45; % °C
T_coolant_in = 30; % °C (realistic for cooling)
T_initial = 25; % °C (room temperature)

%% Cooling System Parameters
% Cooling plate with practical effects
pipe_surface_area = 0.2; % m^2 (multi-pass plate)
h_conv_base = 800; % Base convective coefficient (W/m^2-K)
k_control = 75; % Proportional control gain (W/m^2-K/°C, tuned for practical response)
tau_delay = 30; % Coolant delay time constant (s)

% Practical additions
mass_factor = 1.2; % 20% increase for casing and coolant mass
ambient_noise_std = 0.5; % Ambient temperature noise (±0.5°C)
heat_variation_amplitude = 0.05; % ±5% variation in heat generation
heat_variation_period = 600; % Variation period (s)

%% Simulation
mCp = mass_factor * cell_mass * num_series * num_parallel * cell_specific_heat; % J/K (with inertia)
T_max = NaN; % Initialize
T_final = NaN;
temp_profile = [];
time_data = [];

if use_simulink
    % Create Simulink model
    model_name = 'BatteryCoolingModel';
    bdclose(model_name);
    new_system(model_name);
    % open_system(model_name); % Uncomment for debugging

    % Blocks
    % Integrator: dT/dt -> T
    add_block('simulink/Commonly Used Blocks/Integrator', [model_name '/ThermalDynamics']);
    set_param([model_name '/ThermalDynamics'], ...
        'InitialCondition', num2str(T_initial));

    % Heat generation with variation
    add_block('simulink/Sources/Sine Wave', [model_name '/HeatVariation']);
    set_param([model_name '/HeatVariation'], ...
        'Amplitude', num2str(total_heat_base * heat_variation_amplitude), ...
        'Frequency', num2str(2 * pi / heat_variation_period), ...
        'Phase', '0');

    add_block('simulink/Sources/Constant', [model_name '/HeatBase']);
    set_param([model_name '/HeatBase'], 'Value', num2str(total_heat_base));

    add_block('simulink/Commonly Used Blocks/Sum', [model_name '/SumHeat']);
    set_param([model_name '/SumHeat'], 'Inputs', '++');

    % Ambient temperature with noise
    add_block('simulink/Sources/Random Number', [model_name '/AmbientNoise']);
    set_param([model_name '/AmbientNoise'], ...
        'Mean', '0', ...
        'Variance', num2str(ambient_noise_std^2), ...
        'SampleTime', '1');

    add_block('simulink/Sources/Constant', [model_name '/AmbientBase']);
    set_param([model_name '/AmbientBase'], 'Value', num2str(T_ambient_base));

    add_block('simulink/Commonly Used Blocks/Sum', [model_name '/SumAmbient']);
    set_param([model_name '/SumAmbient'], 'Inputs', '++');

    % Temperature feedback
    add_block('simulink/Commonly Used Blocks/Sum', [model_name '/TempDiff']);
    set_param([model_name '/TempDiff'], 'Inputs', '+-');

    % Proportional cooling gain: k_control * (T - T_ambient)
    add_block('simulink/Commonly Used Blocks/Gain', [model_name '/ControlGain']);
    set_param([model_name '/ControlGain'], 'Gain', num2str(k_control * pipe_surface_area));

    % Coolant delay (first-order low-pass filter)
    add_block('simulink/Continuous/Transfer Fcn', [model_name '/CoolantDelay']);
    set_param([model_name '/CoolantDelay'], ...
        'Numerator', '[1]', ...
        'Denominator', sprintf('[%d 1]', tau_delay));

    % Base cooling: h_conv_base * (T - T_coolant)
    add_block('simulink/Commonly Used Blocks/Gain', [model_name '/BaseCoolingGain']);
    set_param([model_name '/BaseCoolingGain'], 'Gain', num2str(-h_conv_base * pipe_surface_area));

    add_block('simulink/Sources/Constant', [model_name '/CoolantTemp']);
    set_param([model_name '/CoolantTemp'], 'Value', num2str(T_coolant_in));

    add_block('simulink/Commonly Used Blocks/Sum', [model_name '/SumBaseCooling']);
    set_param([model_name '/SumBaseCooling'], 'Inputs', '+-');

    % Nonlinear cooling adjustment
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

    % Combine cooling terms
    add_block('simulink/Commonly Used Blocks/Sum', [model_name '/TotalCooling']);
    set_param([model_name '/TotalCooling'], 'Inputs', '+-');

    % Scale by 1/(m*Cp)
    add_block('simulink/Commonly Used Blocks/Gain', [model_name '/ThermalMass']);
    set_param([model_name '/ThermalMass'], 'Gain', num2str(1/mCp));

    % Scope
    add_block('simulink/Sinks/Scope', [model_name '/TempScope']);
    set_param([model_name '/TempScope'], 'SaveToWorkspace', 'on', ...
        'VariableName', 'TempData');

    % Connections
    add_line(model_name, 'HeatBase/1', 'SumHeat/1', 'autorouting', 'on');
    add_line(model_name, 'HeatVariation/1', 'SumHeat/2', 'autorouting', 'on');
    add_line(model_name, 'AmbientBase/1', 'SumAmbient/1', 'autorouting', 'on');
    add_line(model_name, 'AmbientNoise/1', 'SumAmbient/2', 'autorouting', 'on');
    add_line(model_name, 'ThermalDynamics/1', 'TempDiff/1', 'autorouting', 'on');
    add_line(model_name, 'SumAmbient/1', 'TempDiff/2', 'autorouting', 'on');
    add_line(model_name, 'TempDiff/1', 'ControlGain/1', 'autorouting', 'on');
    add_line(model_name, 'ControlGain/1', 'CoolantDelay/1', 'autorouting', 'on');
    add_line(model_name, 'ThermalDynamics/1', 'BaseCoolingGain/1', 'autorouting', 'on');
    add_line(model_name, 'BaseCoolingGain/1', 'SumBaseCooling/1', 'autorouting', 'on');
    add_line(model_name, 'CoolantTemp/1', 'SumBaseCooling/2', 'autorouting', 'on');
    add_line(model_name, 'ThermalDynamics/1', 'NonlinearGain/1', 'autorouting', 'on');
    add_line(model_name, 'NonlinearGain/1', 'NonlinearFactor/1', 'autorouting', 'on');
    add_line(model_name, 'NonlinearBase/1', 'NonlinearFactor/2', 'autorouting', 'on');
    add_line(model_name, 'NonlinearFactor/1', 'NonlinearCooling/1', 'autorouting', 'on');
    add_line(model_name, 'SumBaseCooling/1', 'ApplyNonlinear/1', 'autorouting', 'on');
    add_line(model_name, 'NonlinearCooling/1', 'ApplyNonlinear/2', 'autorouting', 'on');
    add_line(model_name, 'CoolantDelay/1', 'TotalCooling/1', 'autorouting', 'on');
    add_line(model_name, 'ApplyNonlinear/1', 'TotalCooling/2', 'autorouting', 'on');
    add_line(model_name, 'SumHeat/1', 'ThermalMass/1', 'autorouting', 'on');
    add_line(model_name, 'TotalCooling/1', 'ThermalMass/2', 'autorouting', 'on');
    add_line(model_name, 'ThermalMass/1', 'ThermalDynamics/1', 'autorouting', 'on');
    add_line(model_name, 'ThermalDynamics/1', 'TempScope/1', 'autorouting', 'on');

    % Simulation settings
    set_param(model_name, ...
        'StopTime', '3600', ... % 60 min to observe stabilization
        'Solver', 'ode23t', ... % Stiff solver
        'RelTol', '1e-4', ...
        'AbsTol', '1e-6');

    % Simulate
    try
        sim(model_name);
        temp_var = evalin('base', 'TempData');
        T_max = max(temp_var.signals.values);
        T_final = temp_var.signals.values(end);
        temp_profile = temp_var.signals.values;
        time_data = temp_var.time;
    catch e
        fprintf('Simulation failed: %s\n', e.message);
        T_max = NaN;
        T_final = NaN;
        temp_profile = [];
        time_data = [];
        bdclose(model_name);
    end

    bdclose(model_name);
else
    % Numerical ODE solver
    t = 0:0.01:3600; % 0.01s step for accuracy
    T = zeros(size(t));
    T(1) = T_initial;
    Q_cool_delayed = 0; % Initialize delayed cooling
    dt = 0.01; % Time step
    for j = 2:length(t)
        % Heat generation with variation
        total_heat = total_heat_base * (1 + heat_variation_amplitude * sin(2 * pi * t(j-1) / heat_variation_period));

        % Ambient temperature with noise
        T_ambient = T_ambient_base + ambient_noise_std * randn();

        % Nonlinear cooling: h_conv decreases with temperature difference
        delta_T = T(j-1) - T_coolant_in;
        h_eff = h_conv_base / (1 + 0.01 * delta_T);

        % Proportional cooling: k_control * (T - T_ambient) when T > T_ambient
        h_control = 0;
        if T(j-1) > T_ambient
            h_control = k_control * (T(j-1) - T_ambient);
        end
        h_total = h_eff + h_control;

        % Coolant delay (exponential smoothing)
        Q_cool = h_total * pipe_surface_area * (T(j-1) - T_coolant_in);
        Q_cool_delayed = Q_cool_delayed + (Q_cool - Q_cool_delayed) * (dt / tau_delay);

        % Thermal dynamics
        dT_dt = (total_heat - Q_cool_delayed) / mCp;
        T(j) = T(j-1) + dT_dt * dt;
    end
    T_max = max(T);
    T_final = T(end);
    temp_profile = T;
    time_data = t;
end

%% Results
fprintf('Cell Temperature (Max): %.2f °C\n', T_max);
fprintf('Cell Temperature (Final): %.2f °C\n', T_final);
fprintf('Coolant Temperature: %.2f °C\n', T_coolant_in);
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

%% Plot: Transient Temperature Profile
figure('Name', 'Temperature Profile');
if ~isempty(temp_profile)
    plot(time_data, temp_profile, 'b-', 'LineWidth', 2, 'DisplayName', 'Cell Temp');
    hold on;
end
plot([0 3600], [T_coolant_in T_coolant_in], 'c--', 'LineWidth', 1.5, ...
    'DisplayName', sprintf('Coolant Temp (%.0f°C)', T_coolant_in));
plot([0 3600], [T_ambient_base T_ambient_base], 'm--', 'LineWidth', 1.5, ...
    'DisplayName', sprintf('Ambient Temp (%.0f°C)', T_ambient_base));
plot([0 3600], [45 45], 'r--', 'LineWidth', 1.5, 'DisplayName', '45°C Limit');
xlabel('Time (s)');
ylabel('Temperature (°C)');
title('Cell Temperature Over Time (Practical Model)');
grid on;
legend('Location', 'best');

%% Save Results
save('BatteryCoolingResults.mat', 'T_max', 'T_final', 'T_coolant_in', 'T_ambient_base', ...
    'temp_profile', 'time_data');
