function [entropy_open, entropy_closed, entropy_band_open, entropy_band_closed, entropy_averages] = CalcAvgER()

% Calculate average entropy rate per subject

entropy_source_open = load("entropy_source_open.mat").entropy_source_open;
entropy_source_closed = load("entropy_source_closed.mat").entropy_source_closed;

entropy_open = CalcAvgSubjectEntropy(entropy_source_open);
entropy_closed = CalcAvgSubjectEntropy(entropy_source_closed);

% Calculate average entropy rate per subject per band

entropy_source_band_open = load("entropy_source_band_open.mat").entropy_source_band_open;
entropy_source_band_closed = load("entropy_source_band_closed.mat").entropy_source_band_closed;

entropy_band_open = CalcAvgSubjectBandEntropy(entropy_source_band_open);
entropy_band_closed = CalcAvgSubjectBandEntropy(entropy_source_band_closed);

% Calculate average entropy rate overall and per band

entropy_averages = struct(...
    'entropy_open_average', mean(cell2mat(struct2cell(entropy_open)), "omitmissing"), ...
    'entropy_closed_average', mean(cell2mat(struct2cell(entropy_closed)), "omitmissing"), ...
    'entropy_band_open_average', mean(cell2mat(struct2cell(entropy_band_open)), "omitnan"), ...
    'entropy_band_closed_average', mean(cell2mat(struct2cell(entropy_band_closed)), "omitnan"));

% Helper functions

function [entropy_subject_list] = CalcAvgSubjectEntropy(results)
    fields = fieldnames(results);
    entropy_subject_list = struct();
    
    for i = 1:numel(fields)
        subject = fields{i};
        data = results.(subject);
        average_entropy = mean(data, "omitnan");
        entropy_subject_list.(subject) = average_entropy;
    end
end

function [entropy_subject_band_list] = CalcAvgSubjectBandEntropy(results)
    % Each field is a 5x60 matrix (bands x channels)
    fields = fieldnames(results);
    entropy_subject_band_list = struct();
    
    for i = 1:numel(fields)
        subject = fields{i};
        data = results.(subject); % 5x60
        average_band_entropy = mean(data, 2, "omitnan"); % mean across channels, for each band (5x1)
        entropy_subject_band_list.(subject) = average_band_entropy'; % store as 1x5 row vector
    end
end

% Save results to .mat files

save(fullfile('entropy_open.mat'), 'entropy_open');
save(fullfile('entropy_closed.mat'), 'entropy_closed');
save(fullfile('entropy_band_open.mat'), 'entropy_band_open');
save(fullfile('entropy_band_closed.mat'), 'entropy_band_closed');
fprintf('\n*** FILES SAVED *** \n');

% List entropy_source_open entries which have NaN values
% nan_keys_open = fieldnames(entropy_source_open);
% nan_keys_closed = fieldnames(entropy_source_closed);
% nan_keys_open = nan_keys_open(~cellfun(@(x) all(~isnan(entropy_source_open.(x))), nan_keys_open));
% nan_keys_closed = nan_keys_closed(~cellfun(@(x) all(~isnan(entropy_source_closed.(x))), nan_keys_closed));
% fprintf("\n");
% disp(['OPEN epochs subjects with NaN entries: ', nan_keys_open{:}]);
% disp(['CLOSED epochs subjects with NaN entries: ', nan_keys_closed{:}]);

end