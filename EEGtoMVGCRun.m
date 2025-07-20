function [] = EEGtoMVGCRun()

    % Luka Dubravica, 2025
    % Supervised by Hardik Rajpal and Alberto Liardi

    mvgc_formats = {'brain_source', 'region', 'pca', 'source'};

    for i = 1:length(mvgc_formats)
        mvgc_format = mvgc_formats{i};
        fprintf('\nRunning MVGC computation for format: %s\n', mvgc_format);
        
        % Set options and run the computation
        EEGtoMVGCOptions(false, mvgc_format, true);
        
        % Display completion message
        fprintf('MVGC computation for format "%s" completed.\n', mvgc_format);
    end

end