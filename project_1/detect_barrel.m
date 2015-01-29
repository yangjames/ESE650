function [x, y, d] = detect_barrel(rgbim)
load('model.mat');

% get image
G = fspecial('gaussian',[10 10],10);
rgbim = imfilter(rgbim,G,'same');
im = rgb2ycbcr(rgbim);
[r,c,~] = size(im);

% get pixel values
Y = reshape(im(:,:,1),r*c,1);
Cb = reshape(im(:,:,2),r*c,1);
Cr = reshape(im(:,:,3),r*c,1);
colors = double([Y Cb Cr]);

S = [];
iteration = 1;
while isempty(S) && iteration ~= length(models)
    % calculate probability density of each pixel
    P = zeros(r*c,models{iteration}.model1.num_clusters+models{iteration}.model2.num_clusters);
    DM = zeros(r*c,models{iteration}.model1.num_clusters+models{iteration}.model2.num_clusters);
    for j = 1:models{iteration}.model1.num_clusters
        P(:,j) = compute_gaussian_density(colors,models{iteration}.model1.mean(j,:),models{iteration}.model1.cov{j});
        mean_centered = bsxfun(@minus,colors,models{iteration}.model1.mean(j,:));
        DM(:,j) = sqrt(sum((mean_centered/models{iteration}.model1.cov{j}).*mean_centered,2));
    end
    for j = models{iteration}.model1.num_clusters+1:models{iteration}.model1.num_clusters+models{iteration}.model2.num_clusters
        P(:,j) = compute_gaussian_density(colors,models{iteration}.model2.mean(j-models{iteration}.model1.num_clusters,:),models{iteration}.model2.cov{j-models{iteration}.model1.num_clusters});
        mean_centered = bsxfun(@minus,colors,models{iteration}.model2.mean(j-models{iteration}.model1.num_clusters,:));
        DM(:,j) = sqrt(sum((mean_centered/models{iteration}.model2.cov{j-models{iteration}.model1.num_clusters}).*mean_centered,2));
    end

    [~,idx] = max(DM,[],2);
    size(idx)
    mask = reshape(idx > models{iteration}.model1.num_clusters,r,c);
    se_o = strel('square',25);
    se_c = strel('square',8);
    mask = imclose(mask,se_c);
    mask = imopen(mask,se_o);

    CC = bwconncomp(mask);
    S = regionprops(CC);
    indicator = false(length(S),1);
    BB_size = zeros(length(S),2);
    Area = zeros(length(S),1);
    for j = 1:length(S)
        BB_size(j,:) = S(j).BoundingBox(3:4);
        Area(j) = S(j).Area;
    end
    a_filt = Area >= 2280;
    S = S(a_filt);
    depth = zeros(length(S),1);

    for j = 1:length(S)
        [w,w_idx] = min(S(j).BoundingBox(3:4));
        [h,~] = max(S(j).BoundingBox(3:4));
        if abs(w/h-40/57)/(40/57) < 0.3
            indicator(j) = true;
        end
        if w_idx == 1
            depth(j) = 40*11.35/w;
        else
            depth(j) = 57*11.35/h;
        end
    end
    S = S(indicator);
    depth = depth(indicator);
    
    x = zeros(length(S),1);
    y = zeros(length(S),1);
    d = zeros(length(S),1);
    for j = 1:length(S)
        x(j) = S(j).Centroid(1);
        y(j) = S(j).Centroid(2);
        d(j) = depth(j);
    end
    iteration = iteration+1;
end