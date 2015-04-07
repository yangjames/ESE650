%% clear and close everything
clear all
close all

%% add dijkstra code to path
addpath('matlab')

%% load training data
load('sample_car_paths_1.mat');

%% load aerial map and get dimensions
scale = 1/4;
map_rgb = imresize(imread('aerial_color.jpg'),scale);
im_dims = size(map_rgb);
num_pixels = im_dims(1)*im_dims(2);

%% generate features
fprintf('Generating features...\n')
num_features = 6;
features = zeros(num_features,num_pixels);
features(1,:) = double(reshape(map_rgb(:,:,1)/256,1,num_pixels)); % red channel of raw image
features(2,:) = double(reshape(map_rgb(:,:,2)/256,1,num_pixels)); % green channel of raw image
features(3,:) = double(reshape(map_rgb(:,:,3)/256,1,num_pixels)); % blue channel of raw image
features(4,:) = double(reshape(rgb2gray(map_rgb)/256,1,num_pixels)); % grayscale version of raw image
features(5,:) = double(reshape(edge(rgb2gray(map_rgb),'Canny'),1,num_pixels)); % edges of raw image
features(6,:) = double(reshape(medfilt2(rgb2gray(map_rgb)/256),1,num_pixels)); % high pass filter

%% initialize feature weights
weights = rand(1,num_features);
weights = weights/sum(weights);

%% intialize cost map
fprintf('Creating weighted cost map...\n')
cost_map = 0.1 + exp(reshape(weights*features,im_dims(1),im_dims(2)));

%% initial plot of cost map and paths
figure(2)
clf
cost_plot = imshow(cost_map,[min(min(cost_map)) max(max(cost_map))]);
hold on
actual_path = cell(1,length(paths));
optimal_path = cell(1,length(paths));
for i = 1:length(paths)
    actual_path{i} = plot(round(scale*paths{i}(:,1)),round(scale*paths{i}(:,2)),'g-');
    optimal_path{i} = plot(0,0,'r-');
end

%% start gradient descent
fprintf('Starting gradient descent...\n')
converged = false;
optimal_cost = zeros(1,length(paths));
actual_cost = zeros(1,length(paths));
eta = 0.01;
while ~converged
    for i = 1:length(paths)
        fprintf('Calculating cost for path %d of %d...\n',i, length(paths))
        %% scale the path according to map scale
        path_scaled = unique(round(scale*paths{i}),'rows','stable');
        
        %% get segmented map based on min and max of desired path
        map_dims = [min(path_scaled); max(path_scaled)];
        path = bsxfun(@minus,path_scaled,map_dims(1,:))+1;
        cost_map_segment = cost_map(map_dims(1,2):map_dims(2,2),map_dims(1,1):map_dims(2,1));
                
        %% generate optimal path and calculate cost
        start = path(1,:);
        goal = path(end,:);
        fprintf('Calculating cost to go...\n')
        ctg = dijkstra_matrix(cost_map_segment,goal(2),goal(1));
        fprintf('Running Dijkstra...\n')
        [ip1, jp1] = dijkstra_path(ctg, cost_map_segment, start(2), start(1));
        optimal_cost(i) = sum(cost_map_segment(sub2ind(size(cost_map_segment),ip1,jp1)));
        actual_cost(i) = sum(cost_map(sub2ind(im_dims(1:2),path(:,2),path(:,1))));
        
        %% draw the optimal path
        set(optimal_path{i},'xdata',jp1+map_dims(1,1),'ydata',ip1+map_dims(1,2))
        drawnow
    end
    
    %% TODO: adjust weights
    converged = true;
end
