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

im_L_p = imread([img2Path images2(1).name]);
im_R_p = imread([img3Path images3(1).name]);
matcherMex('push',rgb2gray(im_L_p),rgb2gray(im_R_p));

im_L   = imread([img2Path images2(2).name]);
im_R   = imread([img3Path images3(2).name]);
matcherMex('push',rgb2gray(im_L),rgb2gray(im_R));

% match features
matcherMex('match',2);
matched_features = matcherMex('get_matches',2);
for i = 2:num_images
    im = imread([img2Path images2(i).name]);
    set(im_plot,'CData',im)
    drawnow
end