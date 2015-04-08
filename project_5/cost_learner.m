%% clear and close everything
clear all
close all

%% add dijkstra code to path
addpath('matlab')

%% load training data
load('sample_walking_paths_1.mat');

%% load aerial map and get dimensions
scale = 1/4;
map_rgb = imresize(imread('aerial_color.jpg'),scale);
im_dims = size(map_rgb);
num_pixels = im_dims(1)*im_dims(2);

%% generate features
fprintf('Generating features...\n')
map_ycbcr = rgb2ycbcr(map_rgb);

load('pixel_clusters_lab.mat');
cluster_masks = zeros(nColors,num_pixels);
scaled_kmeans = imresize(pixel_labels,scale);
for i = 1:nColors
    cluster_masks(i,:) = reshape(scaled_kmeans == i,1,num_pixels);
end

num_features = 10+nColors;
features = zeros(num_features,num_pixels);

im_scale = 256;
features(1,:) = double(reshape(map_rgb(:,:,1),1,num_pixels))/im_scale; % red channel of raw image
features(2,:) = double(reshape(map_rgb(:,:,2),1,num_pixels))/im_scale; % green channel of raw image
features(3,:) = double(reshape(map_rgb(:,:,3),1,num_pixels))/im_scale; % blue channel of raw image
features(4,:) = double(reshape(rgb2gray(map_rgb),1,num_pixels))/im_scale; % grayscale version of raw image
features(5,:) = double(reshape(edge(rgb2gray(map_rgb),'Canny'),1,num_pixels)); % edges of raw image
features(6,:) = double(reshape(medfilt2(rgb2gray(map_rgb)),1,num_pixels))/im_scale; % high pass filter
features(7,:) = double(reshape(map_ycbcr(:,:,1),1,num_pixels))/im_scale; % brightness channel of raw image
features(8,:) = double(reshape(map_ycbcr(:,:,2),1,num_pixels))/im_scale; % blue chrominance channel of raw image
features(9,:) = double(reshape(map_ycbcr(:,:,3),1,num_pixels))/im_scale; % red chrominance channel of raw image
features(10,:) = double(reshape(rgb2gray(map_rgb) == rgb2gray(map_rgb(1,1,:)),1,num_pixels));

features(11:end,:) = cluster_masks;

%% initialize feature weights
weights = rand(num_features,1);

%% initial plot of cost map and paths
figure(2)
clf
cost_map = 0.1 + exp(reshape(weights'*features,im_dims(1),im_dims(2)));
min_cost = min(min(cost_map));
max_cost = max(max(cost_map));
cost_plot = imshow((cost_map-min_cost)/(max_cost-min_cost),[0 1]);
hold on
actual_path = cell(1,length(paths));
optimal_path = cell(1,length(paths));
for i = 1:length(paths)
    actual_path{i} = plot(round(scale*paths{i}(:,1)),round(scale*paths{i}(:,2)),'g-');
    optimal_path{i} = plot(0,0,'r-');
end

figure(3)
clf
weight_history_plot = cell(num_features,1);
weight_history = weights;
feature_names = {'red','green','blue','gray','edge','median filter','Y','Cb','Cr'};
weight_colors = rand(num_features,3);
for i = 1:num_features
    weight_history_plot{i} = plot(0,weight_history(i),'-o','Color',weight_colors(i,:));
    hold on
end
%legend(feature_names);
grid on

%% start gradient descent
fprintf('Starting gradient descent...\n')
converged = false;

eta = 0.000001;
iteration = 1;
max_iter = 10;

norm_cost = 0;
norm_cost_prev = Inf;
while ~converged && iteration < max_iter
    % generate cost map
    %fprintf('Creating weighted cost map...\n')
    cost_map = 0.1 + exp(reshape(weights'*features,im_dims(1),im_dims(2)));
    min_cost = min(min(cost_map));
    max_cost = max(max(cost_map));
    set(cost_plot,'cdata',(cost_map-min_cost)/(max_cost-min_cost));
    
    feature_cost_actual = zeros(num_features,1);
    feature_cost_optimal = zeros(num_features,1);
    
    for i = 1:length(paths)
        %% scale the path according to map scale
        path_scaled = unique(round(scale*paths{i}),'rows','stable');
        
        %% get segmented map based on min and max of desired path
        map_dims = [min(path_scaled)-1; max(path_scaled)+1];
        path = bsxfun(@minus,path_scaled,map_dims(1,:))+1;
        cost_map_segment = cost_map(map_dims(1,2):map_dims(2,2),map_dims(1,1):map_dims(2,1));
                
        %% generate optimal path
        start = path(1,:);
        goal = path(end,:);
        %{d
        ctg = dijkstra_matrix(cost_map_segment,goal(2),goal(1));
        [optimal_row, optimal_col] = dijkstra_path(ctg, cost_map_segment, start(2), start(1));
        %}
        %{
        first_step = path(2,:)-path(1,:);
        last_step = path(end,:)-path(end-1,:);
        start = [path(1,:) atan2(first_step(2),first_step(1))];
        goal = [path(end,:) atan2(last_step(2),last_step(1))];
        nav = dijkstra_nonholonomic16(cost_map_segment, goal([2 1 3]), start([2 1 3]));
        [optimal_row, optimal_col, ap, cp] = dijkstra_nonholonomic_path(nav, ...
                                                      start(2), start(1), ...
                                                      start(3), 1.0);
         %}                                         
        %% tally up feature costs
        actual_ind = sub2ind(im_dims(1:2),path_scaled(:,2),path_scaled(:,1));
        optimal_ind = sub2ind(im_dims(1:2),round(optimal_row)+map_dims(1,2)-1,round(optimal_col)+map_dims(1,1)-1);
        
        feature_cost_actual = feature_cost_actual + sum(bsxfun(@times,features(:,actual_ind),cost_map(actual_ind)'),2);
        feature_cost_optimal = feature_cost_optimal + sum(bsxfun(@times,features(:,optimal_ind),cost_map(optimal_ind)'),2);
        %{
        figure(4)
        clf
        min_seg = min(min(cost_map_segment));
        max_seg = max(max(cost_map_segment));
        imshow((cost_map_segment-min_seg)/(max_seg-min_seg),[0 1])
        hold on
        plot(path(:,1),path(:,2),'g-')
        plot(optimal_col,optimal_row,'r-')
        %}
        %% draw the optimal path
        set(optimal_path{i},'xdata',optimal_col+map_dims(1,1)-1,'ydata',optimal_row+map_dims(1,2)-1)
 
    end
    %% adjust weights
    norm_cost = norm(feature_cost_actual - feature_cost_optimal);
    fprintf('delta error: %6.6f\n',abs(norm_cost));
    if abs(norm_cost) < 0.000001
        converged = true;
    else
        weights = weights - eta*(feature_cost_actual - feature_cost_optimal);
                
        weight_history = [weight_history weights];
        norm_cost_prev = norm_cost;
        iteration = iteration + 1;
        for i = 1:num_features
            set(weight_history_plot{i},'xdata',(1:iteration)-1,'ydata',weight_history(i,:));
        end
    end
    
    drawnow
end
save('test_model_2.mat','cost_map','scale')