%% Depression tests
% Set up and load data
depress_files = dir('Depression_Networks/*txt'); % Path to Network files storing N x N ROIs Adjacency Matrices 
T = readtable('sub_details.csv'); %Dataframe that contains subject demographic and diagnosis data
assert(height(T) == length(depress_files));
nb_subjects = length(depress_files);

is_depressed = contains(T.Study, 'Depression'); % Change according to the variable of interest in the dataframe

depress_UI = UI;
depress_UI.design.ui = 1*[is_depressed, ~is_depressed, T.Age, T.Age.^2, contains(T.Sex, 'M')];
depress_UI.exchange.ui = '';
depress_nets = cell([nb_subjects, 1]);
for i=1:nb_subjects
  fname = sprintf('Depression_Networks/%04d_RestEyesClosed_.txt', T.SubjectID(i)); %Edit filenames according to the subject id format
  depress_nets{i} = readmatrix(fname);
end
depress_UI.matrices.ui = depress_nets;

% Test for increases under Depression
depress_UI.contrast.ui = [+1 -1 0 0 0 0];
depress_NBS_pos = CompositeNBStest({depress_UI})

% Test for decreases under Depression
depress_UI.contrast.ui = [-1 +1 0 0 0 0];
depress_NBS_neg = CompositeNBStest({depress_UI})

% Save significant positive clusters, and assert there are no negative clusters
assert(depress_NBS_neg.n == 0);
T = array2table(zeros([0,4]), 'VariableNames', {'From', 'To', 'Tvalue', 'Significant'});
for i=1:5
  for j=1:5
    if i == j continue; end
    T = [T; {ROI{i}, ROI{j}, depress_NBS_pos.test_stat{1}(i,j), depress_NBS_pos.con_mat{1}(i,j)}];
  end
end
writetable(T, 'depress_network_results.csv');

