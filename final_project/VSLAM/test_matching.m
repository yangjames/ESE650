clear all
close all
addpath('../libviso2/matlab/')
%% load data
img2Path = '../dataset/sequences/05/image_2/';
imgType = '*.png';
images2 = dir([img2Path imgType]);

%{
img3Path = '../dataset/sequences/05/image_3/';
images3 = dir([img3Path imgType]);
%}
K2 = [7.188560000000e+02 0.000000000000e+00 6.071928000000e+02;...
    0.000000000000e+00 7.188560000000e+02 1.852157000000e+02;...
    0.000000000000e+00 0.000000000000e+00 1.000000000000e+00];
H2 = [[0 -1 0; 0 0 -1; 1 0 0]' [0 0.06 1.65]'];
R_w = [0 -1 0; 0 0 -1; 1 0 0]';

%{
K3 = [ 7.188560000000e+02 0.000000000000e+00 6.071928000000e+02;...
    0.000000000000e+00 7.188560000000e+02 1.852157000000e+02;...
    0.000000000000e+00 0.000000000000e+00 1.000000000000e+00];
H3 = [eye(3) [0 0.54 1.65]'];
%}
% matching parameters
param.nms_n                  = 2;   % non-max-suppression: min. distance between maxima (in pixels)
param.nms_tau                = 50;  % non-max-suppression: interest point peakiness threshold
param.match_binsize          = 50;  % matching bin width/height (affects efficiency only)
param.match_radius           = 200; % matching radius (du/dv in pixels)
param.match_disp_tolerance   = 1;   % du tolerance for stereo matches (in pixels)
param.outlier_disp_tolerance = 5;   % outlier removal: disparity tolerance (in pixels)
param.outlier_flow_tolerance = 5;   % outlier removal: flow tolerance (in pixels)
param.multi_stage            = 1;   % 0=disabled,1=multistage matching (denser and faster)
param.half_resolution        = 1;   % 0=disabled,1=match at half resolution, refine at full resolution
param.refinement             = 2;   % refinement (0=none,1=pixel,2=subpixel)

% init matcher
matcherMex('init',param);

%% run VSLAM
R = H2(1:3,1:3);
t = H2(:,end);

ground_bounds = [408 185 818 370]; % [xmin ymin xmax ymax]

% initialize feature tracker
im2_1 = imread([img2Path images2(1).name]);
matcherMex('push',rgb2gray(im2_1));

im2_2 = imread([img2Path images2(2).name]);
matcherMex('push',rgb2gray(im2_2));
matcherMex('match',0);
features_12 = matcherMex('get_matches',0)';
matched_indices_12 = matcherMex('get_indices',0)';

im2_3 = imread([img2Path images2(3).name]);
matcherMex('push',rgb2gray(im2_3));
matcherMex('match',0);
features_23 = matcherMex('get_matches',0)';
matched_indices_23 = matcherMex('get_indices',0)';

figure(1)
imshow(im2_2);
hold on
[~,idx2,idx3] = intersect(matched_indices_12(:,2), matched_indices_23(:,1));
plot(features_12(idx2,3), features_12(idx2,4),'r*')
plot(features_23(idx3,1), features_23(idx3,2),'b*')