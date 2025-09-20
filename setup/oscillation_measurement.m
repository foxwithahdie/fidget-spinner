function oscillation_measurement()
    % filename: a string of the filename of the video you want to process
    % the the video is not in the same folder, 
    % fname should include the absolute path to the file.
    filename = 'video/fidgetspinner2.mp4';


    % show_image: a boolean (0 or 1) that determines whether or not the
    % video of the fidget spinner is displayed during processing.
    show_image = 0;

    [t_list, freq_list] = process_fidget_spinner_video(filename, show_image);
    
    number_of_spokes = 3;
    
    freq_list = freq_list / number_of_spokes;

    filename_save = 'proj_resources/fidget_spinner.mat';
    save(filename_save, "t_list", "freq_list");
end