function angular_velocity_measurement()
    %T_window: width of FFT window (in seconds) .7 is usually good
    T_window = .7;
    %q: a parameter used to filter which portions of frequency curve to use
    % (filter is by width of flat portions of curve)
    % value should be in interval [0,1]
    % values closer to 1 are most "exclusive".
    % .6 is a pretty good value to choose
    q = .6;
    %showAnalysis: boolean that turns visualization on/off
    showAnalysis = 0;
    %uses FFT and some filtering tricks to extract angular velocity of fidget
    %spinner as a function of time
    %OUTPUTS:
    %tlist: list of times for measured values of angular frequency
    %omega_list: measured value of angular frequency
    
    video_data = load("proj_resources/video_data.mat");

    [t_list,freq_list] = fidget_spinner_FFT(video_data.avg_pixel_values, video_data.frame_rate, T_window, q, showAnalysis);

    file_name_save = 'proj_resources/fast_fourier_transform_data.mat';
    save(file_name_save, "t_list", "freq_list");
end