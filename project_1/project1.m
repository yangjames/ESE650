clear all
clc

%% import rgb image
a = dir;
file_names = {a.name};
valid_files = ~cellfun(@isempty,regexp(file_names,'^\d+\.\d+\.png\s*$'));
file_names = file_names(valid_files);
for i = 1:length(file_names)
    rgb_im = imread(file_names{i});
    figure(1)
    clf
    imshow(rgb_im)
    title(['File ' num2str(i) ' of ' num2str(length(file_names))])
    %keyboard
    %pause
    coords = find(roipoly);
end
%{
rgb_im = imread('9.0.png');
figure(1)
clf
imshow(rgb_im)
a = roipoly;
%}
%{
%% change color space to LAB
%cform = makecform('srgb2lab');
t_im = rgb2ycbcr(rgb_im);%applycform(rgb_im,cform);

%% reshape the LAB image to color vectors
t_im_mat = [reshape(t_im(:,:,1),1,size(t_im,1)*size(t_im,2));
    reshape(t_im(:,:,2),1,size(t_im,1)*size(t_im,2));
    reshape(t_im(:,:,3),1,size(t_im,1)*size(t_im,2))]';

%% extract CbCr channels
im_cbcr = double(t_im(:,:,2:3));
nrows = size(im_cbcr,1);
ncols = size(im_cbcr,2);
im_cbcr = reshape(im_cbcr,nrows*ncols,2);

nColors = 3;
[cluster_idx, cluster_center] = kmeans(im_cbcr, nColors,'Replicates',3);

%% view each cluster
pixel_labels = reshape(cluster_idx,nrows,ncols);
figure(1)
clf
imshow(pixel_labels,[]),title('image labeled by cluster index');
a = roipoly
%}
%{
%% plot stuff
figure(1)
clf

%[indices,centroids]=kmeans(prin_score,26,'MaxIter',500);

plot(im_mat(:,2),im_mat(:,3),'.')
xlim([0 255])
ylim([0 255])
zlim([0 255])
xlabel('red')
ylabel('green')
zlabel('blue')
grid on
size(im_mat)
%imshow(im)
%}