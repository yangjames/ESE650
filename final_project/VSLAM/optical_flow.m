%clear all
close all
addpath('../libviso2/matlab/')
%% load data
img2Path = '../dataset/sequences/05/image_2/';
imgType = '*.png';
images2 = dir([img2Path imgType]);

K = [7.188560000000e+02 0.000000000000e+00 6.071928000000e+02;...
    0.000000000000e+00 7.188560000000e+02 1.852157000000e+02;...
    0.000000000000e+00 0.000000000000e+00 1.000000000000e+00];
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

ground_bounds = [408 185 818 370]; % [xmin ymin xmax ymax]

%% get seed 3d points
% get matches between first two frames
im2_prev = imread([img2Path images2(1).name]);
im2_rgb = imread([img2Path images2(4).name]);
imsize = size(im2_rgb);
matcherMex('push',rgb2gray(im2_prev));
matcherMex('push',rgb2gray(im2_rgb));

matcherMex('match',0);
features = matcherMex('get_matches',0)';
matched_indices = matcherMex('get_indices',0)';

Mx = features(:,[1 3]);
My = features(:,[2 4]);

mask = all(bsxfun(@gt, features, repmat(ground_bounds(1:2),1,2)) & bsxfun(@lt, features, repmat(ground_bounds(3:4),1,2)),2);
% get rid of outliers
V = 0;
[~, ~, inliers] = GetInliersRANSAC([Mx(:,1) My(:,1)], [Mx(:,2),My(:,2)],0.005,10000);
x1 = [Mx(inliers' & mask,1) My(inliers' & mask,1)];
x2 = [Mx(inliers' & mask,2) My(inliers' & mask,2)];


%%
figure(2)
clf
showMatchedFeatures(im2_prev,im2_rgb,x1,x2)