function measured_vs_analytical()
    fft_data = load("proj_resources/fidget_spinner.mat");
    
    freq_list = fft_data.freq_list;
    t_list    = fft_data.t_list;

    initial_condition = frequency_plot("", 0);
    [a, b, c] = gov_eq_comparison(0);

    [t, sol] = ode45(@(t, freq) (a * freq.^2 + b * freq + c), [0, 40], initial_condition);
    
    frequency_plot("Measured Data vs Analytical Solution", 1); hold on
        plot(t, sol);
        ylim([0, inf]);
        legend("Measured Data", "Analytical Solution");
    hold off
end