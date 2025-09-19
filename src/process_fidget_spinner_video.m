% Extracts angular velocity of fidget spinner from the video.
% INPUTS:
%   filename: a string of the filename of the video you want to process
%           the the video is not in the same folder, fname should include
%           the absolute path to the file
%   optional argument: you can also input whether or not to you want to
%           visualize the FFT (FFT will be shown by default)
%           0: turn off visualization, 1: keep visualization on
% OUTPUTS:
%   t_list: list of times for measured values of angular frequency
%   freq_list: measured value of angular frequency
function [t_list,freq_list] = process_fidget_spinner_video(filename, varargin)
    T_window = 1;
    q = .6;
    showAnalysis = 1;

    if nargin == 2 && varargin{1} == 0
        showAnalysis = 0;
    end

    [data_mat, frame_rate] = extract_data_from_video(filename);
    avg_pixel_values = optimize_coefficients(data_mat);

    [t_list, freq_list] = fidget_spinner_FFT_complex(avg_pixel_values, frame_rate, T_window, q, showAnalysis);
end


function avg_pixel_values = optimize_coefficients(data_mat)
    % Compute bound on number of streams to use based on size of data
    n_points = min([ floorDiv(1e6, size(data_mat, 1)), size(data_mat, 2), 1000 ]);

    data_mat = data_mat - repmat(mean(data_mat), size(data_mat, 1), 1);
    variance_list = var(data_mat);

    var_filter_param = .7;
    max_iter = 200; E_old = 0; E = 1e5; dE_tol = 1e-14;

    [variance_list,max_var_ind] = sort(variance_list,'descend');

    n_points = min(max(cumsum((variance_list-var_filter_param*variance_list(1))>0)),n_points);

    data_mat = data_mat(:,max_var_ind(1:n_points));

    psi_vector = 2*pi*rand(size(data_mat,1),1);

    Mt = [cos(psi_vector),sin(psi_vector)];
    Mb = (Mt'*Mt)\(Mt'*data_mat);
    
    % Used to determine when to print processing update
    prev_percent_complete = -1;

    % Optimize for coefficients iteratively using linear regression
    for n = 1:max_iter

        % Compute progress in percent
        percent_complete = floor(100*n / max_iter);

        % If the progress has incremented by 1%...
        if percent_complete ~= prev_percent_complete
            % then print the current progress
            disp([ 'Coefficient optimization complete: ', num2str(percent_complete), '%' ]);
            prev_percent_complete = percent_complete;
        end

        if abs(E-E_old)<dE_tol
            break
        end
        E_old = E;
        
        Mt = ((Mb * Mb') \ (Mb * data_mat'))';
        Mt = Mt ./ repmat(sqrt(Mt(:, 1).^2 + Mt(:, 2).^2), 1, 2);
        Mb = (Mt' * Mt) \ (Mt' * data_mat);

        E = norm(Mt * Mb - data_mat);
    end
    
    Mt_norm = sqrt(Mt(:, 1).^2 + Mt(:, 2).^2);
    avg_pixel_values = Mt(:,1) ./ Mt_norm + sqrt(-1) * Mt(:, 2) ./ Mt_norm;

    disp('Coefficient optimization is 100% complete');

end


% Uses FFT and some filtering tricks to extract angular velocity of fidget
% spinner as a function of time.
% INPUTS:
%   avg_pixel_values: list of averaged pixel values from video
%   frame_rate: frame rate
%   T_window: width of FFT window (in seconds) .7 is usually good
%   q: a parameter used to filter which portions of frequency curve to use
%       (filter is by width of flat portions of curve)
%       value should be in interval [0,1]
%       values closer to 1 are most "exclusive". 
%       .6 is a pretty good value to choose
%   showAnalysis: boolean that turns visualization on/off
% OUTPUTS:
%   t_list: list of times for measured values of angular frequency
%   omega_list: measured value of angular frequency
function [t_list, omega_list] = fidget_spinner_FFT_complex(avg_pixel_values, frame_rate, T_window, q, showAnalysis)
    % Compute number of data points to use in FFT
    N_window = ceil(frame_rate * T_window);
    
    % Perform windowed FFT of data and pull out peak frequencies
    [t_list1, freq_list1] = extract_max_integer_freq_complex(avg_pixel_values, frame_rate, N_window, showAnalysis);
    
    % Extract the midpoints and heights of the large flat portions of signal
    [midpoint_list, height_list, width_list] = compute_contiguous_centroids(t_list1, freq_list1);
    [t_list, heights_filtered] = filter_by_widths(midpoint_list, height_list, width_list, q);
    
    % De-alias signal using the fact that frequency is decreasing
    omega_list = de_alias_signal(heights_filtered, N_window);

    % %convert from integer frequencies to angular velocity
    omega_factor = 2*pi*frame_rate/N_window;
    omega_list = omega_factor*omega_list;

    % If showAnalysis is true, generate plots showing different stages of
    % the filtered data.
    if showAnalysis
        figure();
        hold on
            plot(t_list1, omega_factor * freq_list1, Color='k', LineWidth=1, DisplayName="Initial Data");
            plot(t_list, omega_factor * heights_filtered,'o', Color='b', MarkerFaceColor='b', MarkerSize=3, DisplayName="Plateau Centroids");
            plot(t_list,omega_list,'o', Color='r', MarkerFaceColor='r',MarkerSize=3, DisplayName="Alias Correction");
            xlabel('Time (s)');
            ylabel('Frequency (rad/sec)');
            title('Filtering Process');
            legend();
        hold off
    end

    t_list = t_list - t_list(1);
end

% Perform windowed FFT on data and exract the maximum integer frequency of
% the FFT for each window.
% INPUTS:
%   avg_pixel_values: list of averaged pixel values from video
%   frame_rate: frame rate
%   N_window: number of data points to use in each windowed FFT
%   showFFT: boolean that turns visualization on/off
% OUTPUTS:
%   t_list: list of times for each FFT
%   freq_list: maximum integer frequency of each FFT
function [t_list, freq_list] = extract_max_integer_freq_complex(avg_pixel_values, frame_rate, N_window, showFFT) 
    % Used for plotting signal being processed
    t_list_base = (0:(length(avg_pixel_values) - 1)) / frame_rate;

    % Shift values of y so it has a mean of 0
    avg_pixel_values = avg_pixel_values - mean(avg_pixel_values);
    % Compute how many points there are in y
    num_samples = length(avg_pixel_values);

    % Initialize the return list
    freq_list = zeros(1, num_samples - N_window + 1);
    
    % Initialize plotting objects as empty arrays
    fft_fig = []; omega_fig = [];
    window_left_plot_1 = []; window_right_plot_1 = [];
    window_left_plot_2 = []; window_right_plot_2 = [];
    omega_plot = []; peak_plot = []; fft_plot = [];

    % If we are plotting, initialize the plots
    if showFFT
        % Create figure showing peak frequencies as a function of time
        omega_fig = figure();
        hold on
        
            omega_plot = plot(0, 0, Color='k', LineWidth=1);
            axis([ 0, t_list_base(end), 0, N_window ]);
            
            title('Peak Frequency of Windowed FFT');
            xlabel('Time (s)');
            ylabel('Peak Integer Frequency (-)');

            % Create figure showing signal and fft of window
            fft_fig=figure(); subplot(3,1,1); hold on
                title('Video Signal');
                xlabel('Time (s)');
                ylabel('real(y(t)) (-)');
                plot(t_list_base, real(avg_pixel_values), Color='k', LineWidth=1);
                axis([ min(t_list_base), max(t_list_base), -1.5, 1.5 ]);
                window_left_plot_1 = plot(0, 0, Color='r', LineWidth=1);
                window_right_plot_1 = plot(0, 0, Color='r',LineWidth=1);
            hold off
    
            subplot(3,1,2); hold on
                title('Video Signal');
                xlabel('Time (s)');
                ylabel('imag(y(t)) (-)');
                plot(t_list_base, imag(avg_pixel_values), Color='k',LineWidth=1);
                axis([ min(t_list_base), max(t_list_base), -1.5, 1.5 ]);
                window_left_plot_2 = plot(0, 0, Color='r', LineWidth=1);
                window_right_plot_2 = plot(0,0,Color='r',LineWidth=1);
            hold off

            subplot(3,1,3); hold on
                title('Windowed FFT of Video Signal');
                xlabel('Integer Frequency (-)');
                ylabel('Normalized FFT Mag ');
                axis([ 0, N_window, 0, 1 ]);
                fft_plot = plot(0, 0, Color='k', LineWidth=1);
                peak_plot = plot(0, 0, 'o', Color='r', MarkerFaceColor='r', MarkerSize=4);
            hold off
        hold off
    end

    % Iterate through data
    for n = 0:(num_samples - N_window)
        % Left and right indices of the window
        index_0 = n + 1; index_1 = n + N_window;
        % Data points in window
        y_vals = avg_pixel_values(index_0:index_1);

        % Compute and normalize FFT data (divide by peak height)
        dfty = abs(fft(y_vals));
        dfty = dfty / max(dfty);

        % Find integer frequency of peak
        [~, max_index] = max(dfty);

        % Store integer frequency in output list
        freq_list(index_0) = max_index - 1;

        % If we are displaying the plots, update the plots
        if showFFT && mod(index_0, 2) == 0
            % Update plot showing window and corresponding FFT
            set(0, CurrentFigure=fft_fig);
            set(window_left_plot_1, XData=t_list_base(index_0) * [1, 1], YData=1.5 * [-1, 1]);
            set(window_right_plot_1, XData=t_list_base(index_1) * [1, 1], YData=1.5 * [-1, 1]);
            set(window_left_plot_2, XData=t_list_base(index_0) * [1, 1], YData=1.5 * [-1, 1]);
            set(window_right_plot_2, XData=t_list_base(index_1) * [1, 1], YData=1.5 * [-1, 1]);
            set(fft_plot, XData=0:length(dfty) - 1,YData=dfty);
            set(peak_plot, XData=max_index - 1,YData=1);
            drawnow;

            % Update plot showing peak frequency as function of time
            set(0, CurrentFigure=omega_fig);
            set(omega_plot, XData=t_list_base(1:index_0), YData=freq_list(1:index_0));
            drawnow;
        end
    end
    % Create t_list based on frame rate and length of freq_list
    t_list = (0:(length(freq_list) - 1)) / frame_rate;
end

% Given a data set (x_i,y_i), computes the locations and width of the flat
% portions of the line plot of (x,y).
% INPUTS:
%   x_list: list of x coordinates of data set (assumed to be monotonic)
%   y_list: list of y coordinates of data set
% OUTPUTS:
%   midpoint_list: midpoint of each flat portion of data set
%   height_list: height of each flat portion of data set
%   width_list: width of each flat portion of data set
function [midpoint_list, height_list, width_list] = compute_contiguous_centroids(x_list, y_list)
    % Initialize output lists
    midpoint_list = []; height_list = []; width_list = [];

    % Initialize current index
    count_1 = 1;
    
    % Iterate until we reach end of list
    while count_1 <= length(y_list)
        % count_2 is 1 + value of rightmost index of a contiguous flat region
        % Initialize count_2
        count_2 = count_1;

        % Increment count)2 until we reach end of flat region
        while count_2 <= length(y_list) && y_list(count_2) == y_list(count_1)
            count_2 = count_2 + 1;
        end
        
        % Compute midpoint location and store it
        midpoint_list(end+1) = (x_list(count_1) + x_list(count_2 - 1)) / 2;

        % Compute height and store it
        height_list(end+1) = y_list(count_1);

        % Compute width and store it
        width_list(end+1) = abs(x_list(count_2 - 1) - x_list(count_1));

        % Set index to be leftmost index of next flat region
        count_1 = count_2;
    end
end

% Given the results of compute_contiguous_centroids, rejects any flat
% candidates that are to oscillatory given their location in the data
% (as time progresses, flat portions should get wider).
% INPUTS:
%   midpoint_list: midpoints of each flat portion of data set
%   height_list: heights of each flat portion of data set
%   width_list: width of each flat portion of data set
%   q: a parameter used to filter out flat portions that are too narrow
%       value should be in interval [0, 1]
%       values closer to 1 are most "exclusive". 
%       .6 is a pretty good value to choose
% OUTPUTS:
%   midpoints_filter: list of midpoints that weren't rejected
%   heights_filtered: list of heights that weren't rejected
function [midpoints_filtered, heights_filtered] = filter_by_widths(midpoint_list, height_list, width_list, q)
    % List of max width value seen up to this point
    max_width_list = 0 * width_list;
    %Initialize max_width_list
    max_width_list(1) = width_list(1);

    % Populate max_width_list
    for n = 2:length(width_list)
        max_width_list(n)=max(max_width_list(n - 1), width_list(n));
    end

    % Initialize return lists
    midpoints_filtered = []; heights_filtered = [];

    % Iterate through data and reject and flat regions in data that are not
    % sufficiently wide given placement in data.
    for n = 1:length(width_list)
        if width_list(n) >= q * max_width_list(min(max(n, 4), length(max_width_list)))
            midpoints_filtered(end+1) = midpoint_list(n);
            heights_filtered(end+1) = height_list(n);
        end
    end
end

% Uses the fact that the angular frequency of the fidget spinner is
% monotonically decreasing to back out the correct peak frequency even,
% when the fidget spinner is rotating multiple times per frame.
% INPUTS:
%   peak_int_freqs: list of peak integer frequencies from windowed FFT
%       (already filtered by looking at flat portions)
%   N_fft_samples: number of sample points used in FFT
% OUTPUTS:
%   peak_int_freqs: "adjusted" list of peak integer frequencies
function peak_int_freqs = de_alias_signal(peak_int_freqs, num_FFT_samples)
    % Iterate backwards through each element of peak_int_freqs
    for n = (length(peak_int_freqs) - 1):-1:1
        % Use aliasing formula to compute the lowest integer frequency that
        %  would alias to y(n), but would still have a larger value that
        %  y(n+1) (assuming that y(n+1) has already been corrected.

        temp_1 = peak_int_freqs(n);
        temp_1 = temp_1 + num_FFT_samples * ceil((peak_int_freqs(n + 1) - temp_1) / num_FFT_samples);
        while temp_1 <= peak_int_freqs(n + 1)
            temp_1 = temp_1 + num_FFT_samples;
        end

        temp_2 = -peak_int_freqs(n);
        temp_2 = temp_2 + num_FFT_samples * ceil((peak_int_freqs(n + 1) - temp_2) / num_FFT_samples);
        while temp_2 <= peak_int_freqs(n + 1)
            temp_2 = temp_2 + num_FFT_samples;
        end

        peak_int_freqs(n) = min(temp_1, temp_2);
    end
end

% Converts the video file of a fidget spinner to a matrix of time signals
% by randomly choosing a set of pixel and then measuring their rgb values
% over time.
% INPUTS:
%   filename: a string of the filename of the video you want to process
%           the the video is not in the same folder, fname should include
%           the absolute path to the file
% OUTPUTS:
%   data_mat: a matrix of data points
%       vertical corresponds to frame index
%       horizontal corresponds to spatial + rgb channel
%       sampled channels are randomly chosen
%   frame_rate: framerate of the video
function  [data_mat, frame_rate] = extract_data_from_video(filename)
    % Create a new VideoReader object to process the data
    vid_obj = VideoReader(filename);

    % Extract number of frames, width, and height of video
    num_frames = vid_obj.NumFrames;
    width = vid_obj.Width;
    height = vid_obj.Height;

    % Extract the video length and framerate
    vid_duration = vid_obj.Duration;
    frame_rate = vid_obj.FrameRate;

    % Set the number of streams to use based on the length of the video
    num_rows = min(floorDiv(3e6, vid_duration), 3 * width * height);

    % Choose which 
    sample_indices = randsample(width * height * 3, num_rows)';

    % Initialize the matrix to store the data
    data_mat = zeros(num_frames, num_rows);

    % Used to determine when to print processing update
    prev_percent_complete = -1;

    % Initialize y to a bunch of zeros with the correct length
    y = zeros(1, round(vid_duration / frame_rate));

    for n = 1:num_frames
        % Compute progress in percent
        percent_complete = floor(100*n / num_frames);

        % If the progress has incremented by 1%...
        if percent_complete ~= prev_percent_complete
            % then print the current progress
            disp([ 'Video processing complete: ', num2str(percent_complete), '%' ]);
            prev_percent_complete = percent_complete;
        end

        % Extract out the next frame from the video
        vidFrame = readFrame(vid_obj);
 
        % Update matrix with desired pixel rgb values
        data_mat(n,:) = vidFrame(sample_indices);
    end

    disp('Video processing is 100% complete');
end
