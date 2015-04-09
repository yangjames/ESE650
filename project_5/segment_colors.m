clear all
close all
map_rgb = imread('aerial_color.jpg');
im_dims = size(map_rgb);
num_pixels = im_dims(1)*im_dims(2);

cform = makecform('srgb2lab');
map_lab = applycform(map_rgb,cform);

cbcr_px = reshape(double(map_lab(:,:,2:3)),prod(im_dims(1:2)),2);
nColors = 7;
[cluster_idx, cluster_center] = kmeans(cbcr_px,nColors,'distance','sqEuclidean',...
    'Replicates',3);
pixel_labels = reshape(cluster_idx,im_dims(1),im_dims(2));
save('pixel_clusters_lab_2.mat','pixel_labels','cluster_center','nColors')