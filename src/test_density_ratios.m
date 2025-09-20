% Testing how different density ratios for the fidget spinner would effect
% the ODE solution that we have come up with.
% INPUTS:
%   density_ratios: Vector of the different values that we are using for
%   the density ratios.
% OUTPUTS:
%   fidget_spinner_rest_vals: The times that the fidget spinner comes to
%   rest.
function [fidget_spinner_rest_vals] = test_density_ratios(density_ratios)
    arguments
        density_ratios (4, 1) double
    end
    fft_data = load("proj_resources/fidget_spinner.mat");
    
    freq_list = fft_data.freq_list;
    t_list    = fft_data.t_list;
    initial_condition = frequency_plot("", 0);
    
    [a, b, c] = gov_eq_comparison(0);

    fidget_spinner_rest_vals = zeros(size(density_ratios));

    for i = 1:size(density_ratios, 1)
        diff_eq = @(t, freq) ((1 / density_ratios(i)) * (a * freq.^2 + b * freq + c));
        [t, sol] = ode45(diff_eq, [0, 150], initial_condition);
        index = find(abs(sol) < 1);
        first_index = index(1);
        fidget_spinner_rest_vals(i) = t(first_index);
    end
    
    
    disp("Fidget Spinner Rest Values");
    for i = 1:4
        disp("Density Ratio: " + density_ratios(i) + ", Rest at Time = " + fidget_spinner_rest_vals(i) + " seconds")
    end
end