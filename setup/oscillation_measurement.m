function oscillation_measurement()
    % filename: a string of the filename of the video you want to process
    % the the video is not in the same folder, 
    % fname should include the absolute path to the file.
    filename = 'video/fidgetspinner2.mp4';

    % window_bounds: a MATLAB struct that indicates the boundaries of the 
    % window to use for averaging the pixel value.
    window_bounds = struct();
    window_bounds.top = 300;
    window_bounds.bottom = 400;
    window_bounds.left = 700;
    window_bounds.right = 800;

    % show_image: a boolean (0 or 1) that determines whether or not the
    % video of the fidget spinner is displayed during processing.

    % set show_image to 1 if you are still trying to figure out
    % the boundaries of the window to use
    % set show_image to 0 to process the video faster
    show_image = 0;
    % Converts the video file of a fidget spinner to a time signal by
    % computing the average pixel value in a window for each frame.

    % OUTPUTS:
    % avg_pixel_values: a list of the averaged pixel value in the window
    % frame_rate: framerate of the video
    [avg_pixel_values, frame_rate] = video_to_signal(filename, window_bounds, show_image);
    
    filename_save = 'proj_resources/video_data.mat';
    save(filename_save, 'avg_pixel_values', 'frame_rate');
end