clear all
close all

addpath('../libviso2/matlab/')

%% load data
imgType = '*.png';
sequence  = '../dataset/sequences/05/';

img2Path = [sequence 'image_2/'];
images2 = dir([img2Path imgType]);
num_images = length(images2);

im = imread([img2Path images2(1).name]);
figure(1)
clf
im_plot = imshow(im);
for i = 2:num_images
    im = imread([img2Path images2(i).name]);
    set(im_plot,'CData',im)
    drawnow
end