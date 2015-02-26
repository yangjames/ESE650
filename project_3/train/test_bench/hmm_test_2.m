clear all
close all
clc

%% orientation plotting initialization
plotting = true;
if plotting
    figure(1)
    clf
    x_plot = plot3(0,0,0,'r-');
    hold on
    y_plot = plot3(0,0,0,'g-');
    z_plot = plot3(0,0,0,'b-');
    axis equal
    grid on
    lims = [-1.1 1.1];
    xlim(lims)
    ylim(lims)
    zlim(lims)
end

%% load models
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

%% begin labelling
file_name = 'beat3_01.txt';
fid = fopen(['../' pattern{pattern_num} '/' file_name]);
data = textscan(fid,'%d %f %f %f %f %f %f',...
    'TreatAsEmpty',{'NA','na'},'CommentStyle','#');
fclose(fid);

time = double(cell2mat(data(:,1)))/1000;
acc = cell2mat(data(:,2:4))/9.81;
gyro = cell2mat(data(:,5:7));

%% initialize
% initial state mean and covariance
x_k = [[1 0 0 0]'; zeros(3,1)];
P_k = eye(6);

% process and sensor noise
Q = diag([78.98;78.98;78.98;18.9;18.9;18.9]);
R = diag([0.5;0.5;0.5;0.001;0.001;0.001]);

tic
start = toc;
orientation = zeros(4,length(time));
orientation(:,1) = [1 0 0 0]';

for i = 2:length(time)
    % get dt
    dt = time(i)-time(i-1);
    v_k = [acc(i,:)'; gyro(i,:)'];
    [x_k, P_k] = ukf(x_k,P_k,Q,R,v_k,dt);

    % get orientation
    orientation = x_k(1:4);

    % calculate probabilities
    for j = 1:length(models)
        labels = 
    end
end