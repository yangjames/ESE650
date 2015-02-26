clear all
close all
clc

%% train
pattern = {'beat3','beat4','circle','eight','inf','wave'};
contents = dir('data');
file_names = cell(length(contents),1);
for file_idx = 1:length(contents)
    file_names{file_idx} = contents(file_idx).name;...
end

num_left_out = 1;
test_set = zeros(length(pattern),num_left_out);
train_set = zeros(length(pattern),5-num_left_out);
for pattern_num = 1:length(pattern)
    acc = [];
    orient = [];
    time = [];
    valid_files = find(~cellfun(@isempty,regexp(file_names,['^' pattern{pattern_num} '.+(\.mat)$'])));
    test_set(pattern_num,:) = floor(rand(1,num_left_out)*5)+1;
    temp = 1:5;
    train_set(pattern_num,:) = temp(~ismember(temp,test_set(pattern_num,:)));
    
    for i = 1:length(train_set(pattern_num,:))
        data = load(['data/' file_names{valid_files(train_set(pattern_num,i))}]);
        acc = [acc data.full_data.acceleration];
        orient = [orient data.full_data.orientation];
        time = [time data.full_data.time'];
    end
    
    M = 47;
    [O,C] = kmeans(orient',M,'maxiter',1000);
    N = 8;
    max_iter = 500;
    params = baum_welch(O,N,max_iter);
    params.C = C;
    fprintf('-----------------------------\n')
    save(['models/' pattern{pattern_num} '_params.mat'],'params')
end


%% test
% load models
contents = dir('models');
model_file_names = cell(length(contents),1);
for file_idx = 1:length(contents)
    model_file_names{file_idx} = contents(file_idx).name;
end
valid_files = find(~cellfun(@isempty,regexp(model_file_names,'.+\.mat')));
models = cell(length(valid_files),1);
for i = 1:length(valid_files)
    models{i} = load(['models/' model_file_names{valid_files(i)}]);
end

% begin test
pattern = {'beat3','beat4','circle','eight','inf','wave'};
margins = cell(length(models),1);
for i = 1:length(models)
    margins{i} = zeros(size(models{i}.params.B,1),num_left_out*6);
end
idx = 1;
for pattern_num = 1:length(pattern)
    addpath(pattern{pattern_num})

    contents = dir(['../' pattern{pattern_num}]);
    file_names = cell(length(contents),1);
    for file_idx = 1:length(contents)
        file_names{file_idx} = contents(file_idx).name;
    end
    valid_files = find(~cellfun(@isempty,regexp(file_names,'.+\.txt')));
    
    for file_idx = 1:length(test_set(pattern_num,:))
        file_name = file_names{valid_files(test_set(pattern_num,file_idx))};
        fid = fopen(['../' pattern{pattern_num} '/' file_name]);
        data = textscan(fid,'%d %f %f %f %f %f %f',...
            'TreatAsEmpty',{'NA','na'},'CommentStyle','#');
        fclose(fid);

        time = double(cell2mat(data(:,1)))/1000;
        acc = cell2mat(data(:,2:4));
        gyro = cell2mat(data(:,5:7));
        O = [acc';gyro'];
        
        %% Viterbi algorithm initialization
        phi = cell(length(models),1);
        labels = zeros(length(models),length(time));
        
        % get first observation
        for i = 1:length(models)
            [~,O_1] = min(sqrt(sum(bsxfun(@minus,O(:,1)',models{i}.params.C).^2,2)));
            phi{i} = zeros(size(models{i}.params.B,1),length(time));
            phi{i}(:,1) = log(models{i}.params.Pi)+log(models{i}.params.B(:,O_1));
        end
        
        % get label wrt each model
        for t = 1:length(time) 
            for m = 1:length(models)
                [~,labels(m,t)] = min(sqrt(sum(bsxfun(@minus,O(:,t)',models{m}.params.C).^2,2)));
            end
        end
        
        figure(pattern_num)
        clf
        for m = 1:length(models)
            subplot(length(models),1,m)
            hold on
            plot(time,labels(m,:),'b-')
            grid on
            drawnow
        end
        
        % calculate probabilities
        logP = zeros(length(models),1);
        for m = 1:length(models)
            for t = 2:length(time)
                phi{m}(:,t) = max(bsxfun(@plus,phi{m}(:,t-1),log(models{m}.params.A)))'...
                    + log(models{m}.params.B(:,labels(m,t)));
            end
            margins{m}(:,idx) = phi{m}(:,end);
            idx = idx+1;
            logP(m) = max(phi{m}(:,end));
        end
        
        % find the winner
        [logprob,pat] = max(logP);
        fprintf(['file: ' file_name ' | pattern: ' pattern{pat} ' | logprog: %6.6f\n'], logprob)
    end
end