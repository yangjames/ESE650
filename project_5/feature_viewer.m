%clear all
close all

scale = 1;
map_rgb = imresize(imread('aerial_color.jpg'),scale);
im_dims = size(map_rgb);
num_pixels = im_dims(1)*im_dims(2);

%% generate features
map_ycbcr = double(rgb2ycbcr(map_rgb));

Y = map_ycbcr(:,:,1);
Cb = map_ycbcr(:,:,2);
Cr = map_ycbcr(:,:,3);

Y_cutoff_low = 80;
Y_cutoff_high = 200;
Y_cut = Y.*(Y>Y_cutoff_low & Y < Y_cutoff_high);
%{
figure(1)
clf
imshow(Y_cut/max(max(Y_cut)),[0 1])
%}
Cb_cutoff_low = 20;
Cb_cutoff_high = 250;
Cb_cut = Cb.*(Cb>Cb_cutoff_low & Cb<Cb_cutoff_high);
%{
figure(2)
clf
imshow(Cb_cut/max(max(Cb_cut)),[0 1])
%}
%% kmeans data
load('pixel_clusters_lab_2.mat')
cform = makecform('srgb2lab');
map_lab = applycform(map_rgb,cform);

cbcr_px = reshape(double(map_lab(:,:,2:3)),prod(im_dims(1:2)),2);

kmeans_data = zeros(3,num_pixels);
idx = [2 3 7];
for i = 1:3
    kmeans_data(i,:) = sum(bsxfun(@minus,cbcr_px,cluster_center(idx(i),:)).^2,2)';
end

%{
distances = zeros(nColors,prod(im_dims(1:2)));
figure(3)
clf
kmeans_plot = imshow(zeros(im_dims(1:2)),[0 1]);
for i = 1:nColors
    distances(i,:) = sum(bsxfun(@minus,cbcr_px,cluster_center(i,:)).^2,2)';
    max_val = max(distances(i,:));
    min_val = min(distances(i,:));
    set(kmeans_plot,'cdata',(reshape(distances(i,:),im_dims(1:2))-min_val)/(max_val-min_val))
    fprintf('current index: %d\n',i)
    drawnow
    pause();
end
%}

%% potential cost map
cost_map = reshape(sum([reshape(Y_cut,1,num_pixels); reshape(Cb_cut,1,num_pixels); -kmeans_data]),im_dims(1:2));

figure(4)
clf
min_val = min(min(cost_map));
max_val = max(max(cost_map));
imshow((cost_map-min_val)/(max_val-min_val),[0 1])