% Plots frequency over time of the fidget spinner. Provides the initial
% condition.
% INPUTS:
%   title_text: Title of the plot
%   show_plot: 1 or 0 for showing plot or not.
% OUTPUTS:
%   initial_condition: The initial condition.
function initial_condition = frequency_plot(title_text, show_plot)
    fft_data = load("proj_resources/fidget_spinner.mat");
    
    freq_list = fft_data.freq_list;
    t_list = fft_data.t_list;
    
    if show_plot
        figure();
        plot(t_list, freq_list, ".");
        title(title_text);
        xlabel("Time (s)"); ylabel("Angular Velocity (rad/s)");
    end

    initial_condition = freq_list(1);
end