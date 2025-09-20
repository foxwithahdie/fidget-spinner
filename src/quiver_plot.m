% Quiver plot of our analytical solution to the differential equation.
% INPUTS:
%   simulated_initial_conditions: Different initial conditions for plots.
function quiver_plot(simulated_initial_conditions)
    arguments
        simulated_initial_conditions (3, 1) double
    end

    fft_data = load("proj_resources/fidget_spinner.mat");
    
    freq_list = fft_data.freq_list;
    t_list    = fft_data.t_list;
    initial_condition = frequency_plot("", 0);
    
    [a, b, c] = gov_eq_comparison(0);

    diff_eq = @(t, freq) (a * freq.^2 + b * freq + c);
    
    [t, sol] = ode45(diff_eq, [0, 40], initial_condition);

    [time_mesh, freq_mesh] = meshgrid(linspace(0, max(t_list), size(t, 1) / 4), linspace(0, max(freq_list) + max(simulated_initial_conditions) / 3, size(sol, 1) / 4));
    
    x_dir = ones(size(time_mesh));
    y_dir = diff_eq(time_mesh, freq_mesh);

    [t_sim_1, sol_sim_1] = ode45(diff_eq, [0, 40], simulated_initial_conditions(1));
    [t_sim_2, sol_sim_2] = ode45(diff_eq, [0, 40], simulated_initial_conditions(2));
    [t_sim_3, sol_sim_3] = ode45(diff_eq, [0, 40], simulated_initial_conditions(3));

    figure();
    quiver(time_mesh, freq_mesh, x_dir, y_dir, DisplayName="Vector Field"); hold on
        ylim([0, inf]);
        xlabel("Time (s)"); ylabel("Angular Velocity (rad/s)");
        title("ODE Compared To Multiple Initial Conditions");
        plot(t_sim_1, sol_sim_1, Color=[1, 0.5, 0], DisplayName="" + simulated_initial_conditions(1) + " rad/s");
        plot(t_sim_2, sol_sim_2, Color=[1, 0.2, 1], DisplayName="" + simulated_initial_conditions(2) + " rad/s");
        plot(t_sim_3, sol_sim_3, Color=[0.1, 1, 0.1], DisplayName="" + simulated_initial_conditions(3) + " rad/s");
        legend();
    hold off
end