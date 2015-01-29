clear all
clc

%% load image data
fprintf('Loading data...\n')
load barrel_nobarrel_data.mat

%% convert barrel pixels to YCbCr and store
fprintf('Gathering pixel data...\n')
data1 = [];
data2 = [];
for i = 1:length(file_names)
    % open image and convert to YCbCr
    im = rgb2ycbcr(imread(file_names{i}));
    
    % extract Y, Cb, and Cr channels
    Y = im(:,:,1);
    Cb = im(:,:,2);
    Cr = im(:,:,3);
    
    % get barrel colors
    Y1 = Y(coords{i});
    Cb1 = Cb(coords{i});
    Cr1 = Cr(coords{i});
    
    % get non-barrel colors
    indices = (1:size(im,1)*size(im,2))';
    indices = indices(~ismember(indices,coords{i}));
    Y2 = Y(indices(1:5:end));
    Cb2 = Cb(indices(1:5:end));
    Cr2 = Cr(indices(1:5:end));
    
    clear im
    data1 = [data1; Y1 Cb1 Cr1];
    data2 = [data2; Y2 Cb2 Cr2];
    clear Y Cb Cr Y1 Cb1 Cr1 Y2 Cb2 Cr2
end
data1 = double(data1);
data2 = double(data2);

%% generate gmm models
model1 = gmm_train(data1,4,0.002);
model2 = gmm_train(data2,2,0.002);

%% plot stuff
figure(1)
clf
sparsity = 100;
plot3(data1(1:sparsity:end,1),data1(1:sparsity:end,2),data1(1:sparsity:end,3),'.')
axis equal
lims = [0 255];
grid on
hold on
xlabel('Y')
ylabel('Cb')
zlabel('Cr')
for i = 1:model2.num_clusters
    % get principle components of covariance matrix
    [V,Lambda] = eig(model2.cov{i});
    R = 1;
    prin_comps = R*sqrt(diag(Lambda));
    
    % make an ellipsoid
    [xe,ye,ze] = ellipsoid(0,0,0,prin_comps(1),prin_comps(2),prin_comps(3));
    
    % rotate axis aligned ellipse
    rot_el = kron(V(:,1),xe) + kron(V(:,2),ye) + kron(V(:,3),ze);
    
    % get ellipse surface points
    rot_x = rot_el(1:size(rot_el,2),:)+model2.mean(i,1);
    rot_y = rot_el(size(rot_el,2)+1:2*size(rot_el,2),:)+model2.mean(i,2);
    rot_z = rot_el(2*size(rot_el,2)+1:end,:)+model2.mean(i,3);
    
    % plot ellipse
    surf(rot_x,rot_y,rot_z);
    alpha(0.1)
end
drawnow

%% test classifications
for i = 1:length(file_names)
    % get image
    rgbim = imread(file_names{i});
    G = fspecial('gaussian',[5 5],4);
    rgbim = imfilter(rgbim,G,'same');
    
    im = rgb2ycbcr(rgbim);
    [r,c,~] = size(im);
    % get pixel values
    Y = reshape(im(:,:,1),r*c,1);
    Cb = reshape(im(:,:,2),r*c,1);
    Cr = reshape(im(:,:,3),r*c,1);
    colors = double([Y Cb Cr]);
    
    % calculate probability density and Mahalanobis distance of each pixel
    P = zeros(r*c,model1.num_clusters+model2.num_clusters);
    DM = zeros(r*c,model1.num_clusters+model2.num_clusters);
    for j = 1:model1.num_clusters
        P(:,j) = compute_gaussian_density(colors,model1.mean(j,:),model1.cov{j});
        mean_centered = bsxfun(@minus,colors,model1.mean(j,:));
        DM(:,j) = sqrt(sum((mean_centered/model1.cov{j}).*mean_centered,2));
    end
    for j = model1.num_clusters+1:model1.num_clusters+model2.num_clusters
        P(:,j) = compute_gaussian_density(colors,model2.mean(j-model1.num_clusters,:),model2.cov{j-model1.num_clusters});
        mean_centered = bsxfun(@minus,colors,model2.mean(j-model1.num_clusters,:));
        DM(:,j) = sqrt(sum((mean_centered/model2.cov{j-model1.num_clusters}).*mean_centered,2));
    end
    
    [~,idx] = max(DM,[],2);
    size(idx)
    mask = reshape(idx > model1.num_clusters,r,c);
    se_o = strel('disk',10);
    se_c = strel('disk',5);
    mask_og = mask;
    mask = imopen(mask,se_o);
    mask = imclose(mask,se_c);
    mask = (mask & mask_og) | mask;
    
    CC = bwconncomp(mask);
    S = regionprops(CC);
    indicator = false(length(S),1);
    BB = zeros(length(S),4);
    Area = zeros(length(S),1);
    for j = 1:length(S)
        BB(j,:) = S(j).BoundingBox;
        A(j) = S(j).Area;
    end
    
    [w,w_idx] = min(BB(:,3:4),[],2);
    [h,h_idx] = max(BB(:,3:4),[],2);
    potential = abs(w./h-40/57) < 0.25 & Area >= 2280;
    
    for j = 1:length(S)
        w = min(S(j).BoundingBox(3:4));
        h = max(S(j).BoundingBox(3:4));
        if abs(w/h-40/57) < 0.25 && S(j).Area >= 2280
            abs(w/h-40/57)
            indicator(j) = true;
        end
    end
    S = S(indicator);
    
    % display image
    figure(2)
    clf
    imshow(rgbim)
    drawnow
    
    % display mask
    figure(3)
    clf
    imshow(mask)
    hold on
    for i = 1:length(S)
        rectangle('Position',S(i).BoundingBox,'EdgeColor','g');
    end
    drawnow
    % wait until a key is pressed
    pause
end