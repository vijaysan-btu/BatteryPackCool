
clear; clc; close all;

%% Create GUI
fig = uifigure('Name', 'Battery Cooling Simulation', ...
    'Position', [100, 100, 600, 500], ...
    'Color', [0.95 0.95 0.95]); % Light gray background

% Main grid layout
grid = uigridlayout(fig, [1, 2], ...
    'ColumnWidth', {'1x', '1x'}, ...
    'RowHeight', {'1x'}, ...
    'Padding', [10 10 10 10], ...
    'ColumnSpacing', 10);

% --- Left Panel: Battery and Drive Cycle Parameters ---
panel_left = uipanel(grid, ...
    'Title', 'Battery Configuration', ...
    'BackgroundColor', [1 1 1], ...
    'FontWeight', 'bold', ...
    'FontSize', 12);
grid_left = uigridlayout(panel_left, [6, 2], ...
    'RowHeight', {'fit', 'fit', 'fit', 'fit', 'fit', 'fit'}, ...
    'ColumnWidth', {'1x', '2x'}, ...
    'Padding', [10 10 10 10], ...
    'RowSpacing', 5);

uilabel(grid_left, 'Text', 'Series Cells:', ...
    'FontWeight', 'bold', 'FontSize', 11);
edit_num_series = uieditfield(grid_left, 'numeric', ...
    'Value', 13, 'Limits', [1, Inf], ...
    'RoundFractionalValues', 'on', ...
    'Tooltip', 'Number of cells in series');

uilabel(grid_left, 'Text', 'Parallel Cells:', ...
    'FontWeight', 'bold', 'FontSize', 11);
edit_num_parallel = uieditfield(grid_left, 'numeric', ...
    'Value', 10, 'Limits', [1, Inf], ...
    'RoundFractionalValues', 'on', ...
    'Tooltip', 'Number of cells in parallel');

uilabel(grid_left, 'Text', 'Nominal Voltage (V):', ...
    'FontWeight', 'bold', 'FontSize', 11);
edit_cell_voltage = uieditfield(grid_left, 'numeric', ...
    'Value', 3.6, 'Limits', [0, Inf], ...
    'Tooltip', 'Voltage per cell (V)');

uilabel(grid_left, 'Text', 'Capacity (Ah):', ...
    'FontWeight', 'bold', 'FontSize', 11);
edit_cell_capacity = uieditfield(grid_left, 'numeric', ...
    'Value', 8.55, 'Limits', [0, Inf], ...
    'Tooltip', 'Capacity per cell (Ah)');

uilabel(grid_left, 'Text', 'C-rates (comma-separated):', ...
    'FontWeight', 'bold', 'FontSize', 11);
edit_c_rates = uieditfield(grid_left, 'text', ...
    'Value', '0.5,2,1,3,0.5,2,1,3,0.5,2,1,3,0.5,2,1,3,0.5,2,1,3,0.5,2,1,3,0.5', ...
    'Tooltip', 'Discharge rates (e.g., 0.5,2,1,3)');

% --- Right Panel: Ambient and Coolant ---
panel_right = uipanel(grid, ...
    'Title', 'Ambient & Coolant', ...
    'BackgroundColor', [1 1 1], ...
    'FontWeight', 'bold', ...
    'FontSize', 12);
grid_right = uigridlayout(panel_right, [4, 2], ...
    'RowHeight', {'fit', 'fit', 'fit', 'fit'}, ...
    'ColumnWidth', {'1x', '2x'}, ...
    'Padding', [10 10 10 10], ...
    'RowSpacing', 5);

uilabel(grid_right, 'Text', 'Nominal Ambient Temp (°C):', ...
    'FontWeight', 'bold', 'FontSize', 11);
edit_T_ambient_base = uieditfield(grid_right, 'numeric', ...
    'Value', 35, 'Limits', [0, 100], ...
    'Tooltip', 'Nominal ambient temperature (°C)');

uilabel(grid_right, 'Text', 'Worst-Case Ambient Temp (°C):', ...
    'FontWeight', 'bold', 'FontSize', 11);
edit_T_ambient_worst = uieditfield(grid_right, 'numeric', ...
    'Value', 40, 'Limits', [0, 100], ...
    'Tooltip', 'Worst-case ambient temperature (°C)');

uilabel(grid_right, 'Text', 'Coolant Type:', ...
    'FontWeight', 'bold', 'FontSize', 11);
coolant_dropdown = uidropdown(grid_right, ...
    'Items', {'Water-Glycol 50/50', 'Water', 'Ethylene Glycol'}, ...
    'Value', 'Water-Glycol 50/50', ...
    'Tooltip', 'Select coolant type');

% Run Simulation Button
run_button = uibutton(grid_right, 'Text', 'Run Simulation', ...
    'ButtonPushedFcn', @(src, event) run_simulation(edit_num_series.Value, edit_num_parallel.Value, ...
    edit_cell_voltage.Value, edit_cell_capacity.Value, edit_c_rates.Value, ...
    edit_T_ambient_base.Value, edit_T_ambient_worst.Value, coolant_dropdown.Value), ...
    'FontWeight', 'bold', 'FontSize', 12, ...
    'BackgroundColor', [0.2, 0.5, 0.8], 'FontColor', [1, 1, 1]);

%% Simulation Function
function run_simulation(num_series, num_parallel, cell_voltage, cell_capacity, c_rates_str, ...
        T_ambient_base, T_ambient_worst, coolant_type)
    % Validate inputs
    try
        if any([num_series, num_parallel, cell_voltage, cell_capacity] <= 0)
            error('Battery configuration, voltage, and capacity must be positive.');
        end
        if any([T_ambient_base, T_ambient_worst] <= 0)
            error('Ambient temperatures must be positive.');
        end
        c_rates = str2num(c_rates_str); %#ok<ST2NM>
        if isempty(c_rates) || any(c_rates <= 0)
            error('C-rates must be a comma-separated list of positive numbers.');
        end
    catch e
        % Display error in dialog box
        result_fig = uifigure('Name', 'Simulation Error', ...
            'Position', [200, 200, 400, 200], ...
            'Color', [0.95 0.95 0.95]);
        uitextarea(result_fig, ...
            'Value', {sprintf('Error: %s', e.message)}, ...
            'Position', [20, 20, 360, 140], ...
            'Editable', 'off', ...
            'BackgroundColor', [1 1 1], ...
            'FontSize', 11);
        return;
    end

    % Fixed parameters
    cell_resistance_base = 0.007; % Ohm, 7 mOhm
    cell_mass = 0.07; % kg
    cell_specific_heat = 900; % J/kg-K
    T_initial = 25; % °C
    T_max_limit = 45; % °C
    T_runaway = 60; % °C
    pipe_surface_area = 0.2; % m²
    mass_factor = 1.2;
    ambient_noise_std = 1.0;

    % Set cooling parameters based on coolant type
    switch coolant_type
        case 'Water-Glycol 50/50'
            Cp_coolant = 3500; % J/kg-K
            h_conv_base = 1000; % W/m²-K
            m_dot = 0.02; % kg/s
            T_coolant_in = 24; % °C
            R_thermal = 0.02; % K/W
            tau_delay = 10; % s
            k_P = 120; k_I = 1.5; k_D = 25;
        case 'Water'
            Cp_coolant = 4180;
            h_conv_base = 1200;
            m_dot = 0.025;
            T_coolant_in = 22;
            R_thermal = 0.015;
            tau_delay = 8;
            k_P = 150; k_I = 2; k_D = 30;
        case 'Ethylene Glycol'
            Cp_coolant = 2400;
            h_conv_base = 800;
            m_dot = 0.015;
            T_coolant_in = 26;
            R_thermal = 0.03;
            tau_delay = 12;
            k_P = 100; k_I = 1; k_D = 20;
    end

    % Scale cooling parameters based on pack size and discharge rate
    pack_size_factor = (num_series * num_parallel) / (13 * 10);
    max_c_rate = max(c_rates);
    h_conv_base = h_conv_base * (1 + 0.1 * (pack_size_factor - 1) + 0.2 * (max_c_rate / 3));
    m_dot = m_dot * (1 + 0.15 * (pack_size_factor - 1) + 0.1 * (max_c_rate / 3));
    T_coolant_in = max(20, T_coolant_in - 2 * (pack_size_factor - 1) - 1 * (max_c_rate / 3));
    R_thermal = max(0.01, R_thermal / (1 + 0.1 * (pack_size_factor - 1)));
    tau_delay = max(5, tau_delay / (1 + 0.1 * (max_c_rate / 3)));
    k_P = k_P * (1 + 0.2 * (pack_size_factor - 1) + 0.15 * (max_c_rate / 3));
    k_I = k_I * (1 + 0.2 * (pack_size_factor - 1) + 0.15 * (max_c_rate / 3));
    k_D = k_D * (1 + 0.2 * (pack_size_factor - 1) + 0.15 * (max_c_rate / 3));
    T_coolant_in = min(T_coolant_in, min(T_ambient_base, T_ambient_worst) - 5);

    % Drive cycle time segments
    n_segments = length(c_rates);
    time_segments = linspace(0, 3600, n_segments);

    % Battery pack calculations
    pack_capacity = num_parallel * cell_capacity;
    rng(42);
    cell_resistance = cell_resistance_base * (1 + 0.1 * (2 * rand(num_series, num_parallel) - 1));
    avg_cell_resistance = mean(cell_resistance(:));

    % Heat generation
    current_total = c_rates * pack_capacity;
    current_cell = current_total / num_parallel;
    heat_gen_cell = zeros(1, length(c_rates));
    for i = 1:length(c_rates)
        heat_gen_cell(i) = (current_cell(i)^2) * avg_cell_resistance + 0.5;
    end
    total_heat = heat_gen_cell * num_series * num_parallel;

    % Thermal mass
    mCp = mass_factor * cell_mass * num_series * num_parallel * cell_specific_heat;

    % Simulate Nominal Case
    t = 0:0.01:3600; dt = 0.01;
    T = zeros(size(t)); T(1) = T_initial;
    Q_cool_delayed = 0; integral_error = 0; prev_error = 0;
    T_coolant_out_sum = 0;
    runaway_flag = false;

    for j = 2:length(t)
        total_heat_t = interp1(time_segments, total_heat, t(j-1), 'previous', 'extrap');
        T_ambient = T_ambient_base + ambient_noise_std * randn();
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
        T_coolant_out_sum = T_coolant_out_sum + Q_cool;
        if T(j) > T_runaway
            runaway_flag = true;
        end
    end
    T_max = max(T);
    T_final = T(end);
    temp_profile = T;
    time_data = t;
    T_coolant_out = T_coolant_in + (T_coolant_out_sum / length(t)) / (m_dot * Cp_coolant);

    % Simulate Worst-Case
    T = zeros(size(t)); T(1) = T_initial;
    Q_cool_delayed = 0; integral_error = 0; prev_error = 0;
    runaway_flag_worst = false;
    for j = 2:length(t)
        total_heat_t = interp1(time_segments, total_heat, t(j-1), 'previous', 'extrap');
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

    % Display Results in Dialog Box
    result_fig = uifigure('Name', 'Simulation Results', ...
        'Position', [200, 200, 800, 500], ...
        'Color', [0.95 0.95 0.95]);
    result_grid = uigridlayout(result_fig, [2, 1], ...
        'RowHeight', {'1x', '1x'}, ...
        'ColumnWidth', {'1x'}, ...
        'Padding', [10 10 10 10], ...
        'RowSpacing', 10);

    % Text area for results
    result_text = uitextarea(result_grid, ...
        'Editable', 'off', ...
        'BackgroundColor', [1 1 1], ...
        'FontSize', 11);
    
    % Axes for plot
    result_axes = uiaxes(result_grid, ...
        'BackgroundColor', [1 1 1], ...
        'FontSize', 10);
    result_axes.XLabel.String = 'Time (s)';
    result_axes.YLabel.String = 'Temperature (°C)';
    result_axes.Title.String = 'Cell Temperature (Nominal Case)';
    result_axes.GridLineStyle = '-';
    grid(result_axes, 'on');

    % Format results
    result_str = {};
    result_str{end+1} = sprintf('Nominal Case (Ambient = %.0f°C):', T_ambient_base);
    result_str{end+1} = sprintf('  Cell Temperature (Max): %.2f °C', T_max);
    result_str{end+1} = sprintf('  Cell Temperature (Final): %.2f °C', T_final);
    result_str{end+1} = sprintf('  Coolant Inlet Temperature: %.2f °C', T_coolant_in);
    result_str{end+1} = sprintf('  Coolant Outlet Temperature (Avg): %.2f °C', T_coolant_out);
    if T_max <= T_max_limit && ~isnan(T_max)
        result_str{end+1} = sprintf('  Cell temperature is within limit (%.0f°C).', T_max_limit);
        if abs(T_final - T_ambient_base) <= 1
            result_str{end+1} = sprintf('  Cell temperature stabilized within ambient range (%.0f ± 1°C).', T_ambient_base);
        else
            result_str{end+1} = sprintf('  Cell temperature stabilized outside ambient range (%.0f ± 1°C) but within limit.', T_ambient_base);
        end
    else
        result_str{end+1} = sprintf('  Warning: Cell temperature exceeds limit (%.0f°C) or simulation failed.', T_max_limit);
    end
    if runaway_flag
        result_str{end+1} = sprintf('  Warning: Thermal runaway risk detected: Temperature exceeded %.0f°C.', T_runaway);
    end

    result_str{end+1} = sprintf('\nWorst-Case (Ambient = %.0f°C):', T_ambient_worst);
    result_str{end+1} = sprintf('  Cell Temperature (Max): %.2f °C', T_max_worst);
    result_str{end+1} = sprintf('  Cell Temperature (Final): %.2f °C', T_final_worst);
    result_str{end+1} = sprintf('  Coolant Inlet Temperature: %.2f °C', T_coolant_in);
    if T_max_worst <= T_max_limit && ~isnan(T_max_worst)
        result_str{end+1} = sprintf('  Cell temperature is within limit (%.0f°C).', T_max_limit);
        if abs(T_final_worst - T_ambient_worst) <= 1
            result_str{end+1} = sprintf('  Cell temperature stabilized within ambient range (%.0f ± 1°C).', T_ambient_worst);
        else
            result_str{end+1} = sprintf('  Cell temperature stabilized outside ambient range (%.0f ± 1°C) but within limit.', T_ambient_worst);
        end
    else
        result_str{end+1} = sprintf('  Warning: Cell temperature exceeds limit (%.0f°C) or simulation failed.', T_max_limit);
    end
    if runaway_flag_worst
        result_str{end+1} = sprintf('  Warning: Thermal runaway risk detected: Temperature exceeded %.0f°C.', T_runaway);
    end
    result_text.Value = result_str;

    % Plot Nominal Case
    cla(result_axes);
    plot(result_axes, time_data, temp_profile, 'b-', 'LineWidth', 2, 'DisplayName', 'Cell Temp (Nominal)');
    hold(result_axes, 'on');
    plot(result_axes, [0 3600], [T_coolant_in T_coolant_in], 'c--', 'LineWidth', 1.5, ...
        'DisplayName', sprintf('Coolant Inlet (%.0f°C)', T_coolant_in));
    plot(result_axes, [0 3600], [T_ambient_base T_ambient_base], 'm--', 'LineWidth', 1.5, ...
        'DisplayName', sprintf('Ambient (%.0f°C)', T_ambient_base));
    plot(result_axes, [0 3600], [45 45], 'r--', 'LineWidth', 1.5, 'DisplayName', '45°C Limit');
    plot(result_axes, [0 3600], [60 60], 'k:', 'LineWidth', 1.5, 'DisplayName', 'Runaway Threshold (60°C)');
    legend(result_axes, 'Location', 'best', 'FontSize', 10);
    hold(result_axes, 'off');
end

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

