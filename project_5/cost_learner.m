%% clear and close everything
clear all
close all

%% add dijkstra code to path
addpath('matlab')

%% load training data
load('car_path_training_set_6.mat');

%% load aerial map and get dimensions
scale = 1/4;
map_rgb = imresize(imread('aerial_color.jpg'),scale);
im_dims = size(map_rgb);
num_pixels = im_dims(1)*im_dims(2);

%% generate features
fprintf('Generating features...\n')

%{d
addpath('feature_scripts_2')

[~,f1] = blueish_structures(map_rgb);
[~,f2] = bright_roofs(map_rgb);
[~,f3] = empty_border(map_rgb);
[~,f4] = grass_trees_sidewalks(map_rgb);
[~,f5] = right_grey_roof(map_rgb);
[~,f6] = shadows(map_rgb);
[~,f7] = sidewalk(map_rgb);
[~,f8] = some_structures(map_rgb);
[~,f9] = trees(map_rgb);
[~,f10] = trees_grass(map_rgb);
[~,f11] = white_border(map_rgb);
[~,f12] = buildings(map_rgb);
[~,f13] = dirt_grass(map_rgb);
[~,f14] = top_left_bright_roof(map_rgb);
[~,f15] = med_grey_roofs(map_rgb);
[~,f16] = building_corners(map_rgb);
[~,f17] = light_shadows(map_rgb);
f18 = abs(gradient(im2double(rgb2gray(map_rgb))));

num_features = 52;
features = zeros(num_features,num_pixels);
features(1:3,:) = im2double(reshape(f1,num_pixels,3)');
features(4:6,:) = im2double(reshape(f2,num_pixels,3)');
features(7:9,:) = im2double(reshape(f3,num_pixels,3)');
features(10:12,:) = im2double(reshape(f4,num_pixels,3)');
features(13:15,:) = im2double(reshape(f5,num_pixels,3)');
features(16:18,:) = im2double(reshape(f6,num_pixels,3)');
features(19:21,:) = im2double(reshape(f7,num_pixels,3)');
features(22:24,:) = im2double(reshape(f8,num_pixels,3)');
featuers(25:27,:) = im2double(reshape(f9,num_pixels,3)');
features(28:30,:) = im2double(reshape(f10,num_pixels,3)');
features(31:33,:) = im2double(reshape(f11,num_pixels,3)');
features(34:36,:) = im2double(reshape(f12,num_pixels,3)');
features(37:39,:) = im2double(reshape(f13,num_pixels,3)');
features(40:42,:) = im2double(reshape(f14,num_pixels,3)');
features(43:45,:) = im2double(reshape(f15,num_pixels,3)');
features(46:48,:) = im2double(reshape(f16,num_pixels,3)');
features(49:51,:) = im2double(reshape(f17,num_pixels,3)');
features(52,:) = im2double(reshape(f18,num_pixels,1)');
%}
%{
map_ycbcr = rgb2ycbcr(map_rgb);
map_hsv = rgb2hsv(map_rgb);

gaussian_filter = fspecial('gaussian',3,2);
map_ycbcr_blur = imfilter(map_ycbcr,gaussian_filter);
map_hsv_blur = imfilter(map_hsv,gaussian_filter);
map_rgb_blur = imfilter(map_rgb,gaussian_filter);

%[~,borders] = empty_border(map_rgb);
grad = gradient(double(rgb2gray(map_rgb)));
[~,roofs] = bright_roofs(map_rgb);

num_features = 16;
features = zeros(num_features,num_pixels);
features(1:3,:) = double(reshape(map_rgb,num_pixels,3)')/256;
features(4:6,:) = double(reshape(map_ycbcr,num_pixels,3)')/256;
features(7:9,:) = double(reshape(map_hsv,num_pixels,3)')/256;
%features(10:12,:) = double(reshape(map_rgb_blur,num_pixels,3)')/256;
%features(13:15,:) = double(reshape(map_ycbcr_blur,num_pixels,3)')/256;
%features(16:18,:) = double(reshape(map_hsv_blur,num_pixels,3)')/256;
features(10:12,:) = double(reshape(map_hsv.^2,num_pixels,3)')/256^2;
%features(22:24,:) = double(reshape(borders,num_pixels,3)')/256;
features(13,:) = double(reshape(grad,num_pixels,1)')/256;
features(14:16,:) = double(reshape(roofs,num_pixels,3)')/256;
%}

%% initialize feature weights
%{d
weights = rand(num_features,1);
%}
%weights = ones(num_features,1);
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
drawnow
%% start gradient descent
fprintf('Starting gradient descent...\n')
converged = false;

eta = 1e-9;
iteration = 1;
max_iter = 10000;

norm_cost = 0;
norm_cost_prev = Inf;
while ~converged && iteration <= max_iter
    % generate cost map
    %fprintf('Creating weighted cost map...\n')
    %disp(weights)
    cost_map = 1e-5 + exp(reshape(weights'*features,im_dims(1),im_dims(2)));
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
        
        ctg = dijkstra_matrix(cost_map_segment,goal(2),goal(1));
        [optimal_row, optimal_col] = dijkstra_path(ctg, cost_map_segment, start(2), start(1));
        
        %% tally up feature costs
        actual_ind = sub2ind(im_dims(1:2),path_scaled(:,2),path_scaled(:,1));
        optimal_ind = sub2ind(im_dims(1:2),round(optimal_row)+map_dims(1,2)-1,round(optimal_col)+map_dims(1,1)-1);
        
        feature_cost_actual = feature_cost_actual + sum(bsxfun(@times,features(:,actual_ind),cost_map(actual_ind)'),2);
        feature_cost_optimal = feature_cost_optimal + sum(bsxfun(@times,features(:,optimal_ind),cost_map(optimal_ind)'),2);

        %% draw the optimal path
        set(optimal_path{i},'xdata',optimal_col+map_dims(1,1)-1,'ydata',optimal_row+map_dims(1,2)-1)
 
    end
    %% adjust weights
    norm_cost = norm(feature_cost_actual - feature_cost_optimal);
    fprintf('delta error: %6.6f\n',norm_cost-norm_cost_prev);
    if abs(norm_cost-norm_cost_prev) < 0.001
        converged = true;
        cost_map = 0.1 + exp(reshape(weights'*features,im_dims(1),im_dims(2)));
    else
        weights = weights - eta*(feature_cost_actual - feature_cost_optimal);
        
        %std_weights = std(weights);
        %mean_weights = mean(weights);
        %weights(end) = (weights(end) - mean_weights)/std_weights;
        
        weight_history = [weight_history weights];
        norm_cost_prev = norm_cost;
        iteration = iteration + 1;
        for i = 1:num_features
            set(weight_history_plot{i},'xdata',(1:iteration)-1,'ydata',weight_history(i,:));
        end
    end
    
    drawnow
end
save('car_model_6.mat','cost_map','scale','weights')