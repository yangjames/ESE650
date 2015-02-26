clear all
close all
%clc

%% load models
model_dir = 'models';
best = 'decent_model_4';
model_dir = best;
contents = dir(model_dir);
model_file_names = cell(length(contents),1);
for file_idx = 1:length(contents)
    model_file_names{file_idx} = contents(file_idx).name;
end
valid_files = find(~cellfun(@isempty,regexp(model_file_names,'.+\.mat')));
models = cell(length(valid_files),1);
for i = 1:length(valid_files)
    models{i} = load([model_dir '/' model_file_names{valid_files(i)}]);
end

%% begin test
pattern = {'beat3','beat4','circle','eight','inf','wave'};
%for pattern_num = 1:length(pattern)
    %addpath(pattern{pattern_num})

    %contents = dir(['../' pattern{pattern_num}]);
    contents = dir('project3_test_data/single');
    file_names = cell(length(contents),1);
    for file_idx = 1:length(contents)
        file_names{file_idx} = contents(file_idx).name;
    end
    valid_files = find(~cellfun(@isempty,regexp(file_names,'.+\.txt')));
    
    for file_idx = 1:length(valid_files)
        file_name = file_names{valid_files(file_idx)};
        %fid = fopen(['../' pattern{pattern_num} '/' file_name]);
        fid = fopen(['project3_test_data/single/' file_name]);
        data = textscan(fid,'%d %f %f %f %f %f %f',...
            'TreatAsEmpty',{'NA','na'},'CommentStyle','#');
        fclose(fid);

        time = double(cell2mat(data(:,1)))/1000;
        acc = cell2mat(data(:,2:4));
        gyro = cell2mat(data(:,5:7));
        %{
        %% run UKF
        % initial state mean, covariance, process noise, and sensor noise
        x_k = [[1 0 0 0]'; zeros(3,1)];
        P_k = eye(6);
        Q = diag([78.98;78.98;78.98;18.9;18.9;18.9]);
        R = diag([0.5;0.5;0.5;0.001;0.001;0.001]);
        
        O = zeros(4,length(time));
        O(:,1) = [1 0 0 0]';
        for i = 2:length(time)
            % get dt
            dt = time(i)-time(i-1);
            v_k = [(acc(i,:)/9.81)'; gyro(i,:)'];
            [x_k, P_k] = ukf(x_k,P_k,Q,R,v_k,dt);
            
            % get orientation
            O(:,i) = x_k(1:4);
        end 
        %}
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
        
        %{
        figure(pattern_num)
        clf
        for m = 1:length(models)
            subplot(length(models),1,m)
            hold on
            plot(time,labels(m,:),'b-')
            grid on
            title(pattern{m})
            drawnow
        end
        %}
        
        % calculate probabilities
        logP = zeros(length(models),1);
        for m = 1:length(models)
            for t = 2:length(time)
                phi{m}(:,t) = max(bsxfun(@plus,phi{m}(:,t-1),log(models{m}.params.A)))'...
                    + log(models{m}.params.B(:,labels(m,t)));
            end
            logP(m) = max(phi{m}(:,end));
        end
        [v,p] = sort(logP);
        
        % find the winner
        [logprob,pat] = max(logP);
        fprintf(['file: ' file_name ' | pattern: ' pattern{p(end)} ', ' pattern{p(end-1)} ', ' pattern{p(end-2)} ' | confidence: %6.6f\n'], -(v(end)-v(end-1))/(v(end)+v(end-1))*100)%' | logprog: %6.6f, %6.6f, %6.6f\n'], v(end),v(end-1),v(end-2))
    end
%end