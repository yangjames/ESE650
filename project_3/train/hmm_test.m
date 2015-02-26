clear all
close all
clc

%% get model files
model_dir = dir('models');
file_names = cell(length(model_dir),1);
for file_idx = 1:length(model_dir)
    file_names{file_idx} = model_dir(file_idx).name;...
end

patterns = {'beat3','beat4','circle','eight','inf','wave'};
valid_files = find(~cellfun(@isempty,regexp(file_names,'.+(\.mat)$')));
model = cell(length(valid_files),1);
for model_file = 1:length(valid_files)
    model{model_file} = load(['models/' file_names{valid_files(model_file)}]);
end

%%