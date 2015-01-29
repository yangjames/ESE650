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
    Y2 = Y(indices(1:3:end));
    Cb2 = Cb(indices(1:3:end));
    Cr2 = Cr(indices(1:3:end));
    
    clear im
    data1 = [data1; Y1 Cb1 Cr1];
    data2 = [data2; Y2 Cb2 Cr2];
    clear Y Cb Cr Y1 Cb1 Cr1 Y2 Cb2 Cr2
end
data1 = double(data1);
data2 = double(data2);

%% generate gmm models
model1 = gmm_train(data1,7,0.002);
model2 = gmm_train(data2,2,0.002);

%% plot stuff
figure(1)
clf
sparsity = 100;
axis equal
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

plot3(data2(1:sparsity:end,1),data2(1:sparsity:end,2),data2(1:sparsity:end,3),'.')
drawnow

%% test classifications
for i = 1:length(file_names)
    % get image
    rgbim = imread(file_names{i});
    G = fspecial('gaussian',[10 10],10);
    rgbim = imfilter(rgbim,G,'same');
    
    im = rgb2ycbcr(rgbim);
    [r,c,~] = size(im);
    % get pixel values
    Y = reshape(im(:,:,1),r*c,1);
    Cb = reshape(im(:,:,2),r*c,1);
    Cr = reshape(im(:,:,3),r*c,1);
    colors = double([Y Cb Cr]);
    
    % calculate probability density of each pixel
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
    se_o = strel('square',25);
    se_c = strel('square',8);
    mask_og = mask;
    %mask = imopen(mask,se_o);
    mask = imclose(mask,se_c);
    mask = imopen(mask,se_o);
    %mask = (mask & mask_og) | mask;
    
    CC = bwconncomp(mask);
    S = regionprops(CC);
    indicator = false(length(S),1);
    BB_size = zeros(length(S),2);
    BB_coords = zeros(length(S),4);
    Area = zeros(length(S),1);
    for j = 1:length(S)
        BB_size(j,:) = S(j).BoundingBox(3:4);
        Area(j) = S(j).Area;
    end
    a_filt = Area >= 2280;
    S = S(a_filt);
    BB_size = BB_size(a_filt,:);
    Area = Area(a_filt);
    depth = zeros(length(S),1);
    %{
    [w,w_idx] = min(BB_size,[],2);
    [h,h_idx] = max(BB_size,[],2);
    potential = abs(w./h-40/57) < 0.25 & Area >= 2280;
    %}
    for j = 1:length(S)
        [w,w_idx] = min(S(j).BoundingBox(3:4));
        [h,h_idx] = max(S(j).BoundingBox(3:4));
        if abs(w/h-40/57)/(40/57) < 0.25
            indicator(j) = true;
        end
        if w_idx == 1
            depth(j) = 40*11.35/w;
        else
            depth(j) = 57*11.35/h;
        end
    end
    S = S(indicator);

    % display image
    figure(2)
    clf
    imshow(rgbim)
    hold on
    for j = 1:length(S)
        plot(S(j).Centroid(1),S(j).Centroid(2),'g+')
        rectangle('Position',S(j).BoundingBox,'EdgeColor','g');
        text(S(j).Centroid(1)+20,S(j).Centroid(2)-20,[num2str(depth(j),5) 'm'],'BackgroundColor','w');
    end
    title(['Image ' num2str(i) ' of ' num2str(length(file_names)) ' ' file_names{i}])
    drawnow
    %{
    % display mask
    figure(3)
    clf
    imshow(mask)
    hold on
    for j = 1:length(S)
        rectangle('Position',S(j).BoundingBox,'EdgeColor','g');
    end
    drawnow
    %}
    % wait until a key is pressed
    pause
end