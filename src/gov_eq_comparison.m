% A plot of d(omega(t))/dt compared to omega(t). Should be quadratic
% INPUTS:
%   show_plot: 0 or 1 for showing the plot.
% OUTPUTS:
%   a, b, c: coefficients for the quadratic fitting plot.
function [a, b, c] = gov_eq_comparison(show_plot)
    fft_data = load("proj_resources/fidget_spinner.mat");
    
    freq_list = fft_data.freq_list;
    t_list = fft_data.t_list;
    
    figure();
    
    angular_acceleration = diff(freq_list) ./ diff(t_list);
    midpoint_velocity = (freq_list(2:end) + freq_list(1:end-1)) / 2;
    coeffs = polyfit(midpoint_velocity, angular_acceleration, 2);
    estimated_velocity = linspace(min(midpoint_velocity), max(midpoint_velocity), 100);
    estimated_accel = polyval(coeffs, estimated_velocity);
    
    if show_plot
        plot(midpoint_velocity, angular_acceleration, '.', DisplayName="Angular Acceleration Comparison"); hold on
            plot(estimated_velocity, estimated_accel, DisplayName="Quadratic Fit")
            title("Fitting Frequency Data to a Quadratic Curve");
            xlabel("Angular Velocity (rad/s)"); ylabel("Angular Acceleration (rad/s^2)");
            legend();
        hold off;
    end

    a = coeffs(1); b = coeffs(2); c = coeffs(3);
end