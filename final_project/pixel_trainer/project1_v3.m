clear all
close all
clc

%% load image data
fprintf('Loading data...\n')
load road_colors.mat

file_path = '../dataset/sequences/00/image_2/';
a = dir(file_path);
file_names = {a.name};
valid_files = ~cellfun(@isempty,regexp(file_names,'^\d+\.png\s*$'));
file_names = file_names(valid_files);

training = false;
if training
%% convert barrel pixels to YCbCr and store
fprintf('Gathering pixel data...\n')
data1 = [];
data2 = [];
num_valid = sum(~cellfun(@isempty,colors));
for i = 1:num_valid
    % open image and convert to YCbCr
    im = rgb2hsv(imread([file_path file_names{i}]));
    
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
model1 = gmm_train(data1,1,0.002);
else
    load('model1.mat');
end
%% plot stuff
%{
figure(1)
clf
sparsity = 1;
axis equal
grid on
hold on
xlabel('Y')
ylabel('Cb')
zlabel('Cr')

plot3(data2(1:sparsity:end,1),data2(1:sparsity:end,2),data2(1:sparsity:end,3),'.')
drawnow
%}

%% test classifications
for i = 1:length(file_names)
    % get image
    rgbim = imread([file_path file_names{i}]);
    %G = fspecial('gaussian',[10 10],10);
    %rgbim = imfilter(rgbim,G,'same');
    
    im = rgb2hsv(rgbim);
    [r,c,~] = size(im);
    % get pixel values
    Y = reshape(im(:,:,1),r*c,1);
    Cb = reshape(im(:,:,2),r*c,1);
    Cr = reshape(im(:,:,3),r*c,1);
    colors = double([Y Cb Cr]);
    
    % calculate probability density of each pixel
    P = zeros(r*c,model1.num_clusters);
    p = zeros(r*c,model1.num_clusters);
    DM = zeros(r*c,model1.num_clusters);
    for j = 1:model1.num_clusters
        P(:,j) = compute_gaussian_density(colors,model1.mean(j,:),model1.cov{j});
        %fprintf('calculating cdf...\n')
        %p(:,j) = mvncdf(colors,model1.mean(j,:),model1.cov{j});
        mean_centered = bsxfun(@minus,colors,model1.mean(j,:));
        DM(:,j) = sqrt(sum((mean_centered/model1.cov{j}).*mean_centered,2));
    end
    %[~,idx] = max(DM,[],2);
    idx = any(P>9,2);
    mask = reshape(idx,r,c);
    se_o = strel('square',5);
    se_c = strel('square',5);
    %mask_og = mask;
    %mask = imopen(mask,se_o);
    %mask = imclose(mask,se_c);
    %mask = imopen(mask,se_o);
    %mask = (mask & mask_og) | mask;

    % display image
    figure(2)
    clf
    imshow(rgbim)
    
    figure(3)
    clf
    imshow(mask)
    
    drawnow
    % wait until a key is pressed
    pause
end