% Compares the ratio of rest times to the density ratios for each type of
% material.
% INPUTS:
%   rest_times: The array of different rest times for the different density
%   ratios.
%   density_ratios: The array of different density ratios for the new
%   against the old material.
function rest_vs_density_ratios(rest_times, density_ratios)
    arguments
        rest_times (4, 1) double
        density_ratios (4, 1) double
    end

    rest_time_ratios = rest_times / (rest_times(2));
    plot(density_ratios, rest_time_ratios, Color=[1, 0.1, 0.1]);
    title("Density Ratios Compared to Rest Time Ratios");
    xlabel("Density Ratio (ρ_{new} / ρ_{old})");
    ylabel("Rest Time Ratio (t_{stop} / t^{*}_{stop}");
end