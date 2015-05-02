clear all
close all

addpath('../libviso2/matlab/')

%% load data
imgType = '*.png';
sequence  = '../dataset/sequences/05/';

img2Path = [sequence 'image_2/'];
images2 = dir([img2Path imgType]);

img3Path = [sequence 'image_3/'];
images3 = dir([img3Path imgType]);

% camera 2 and 3 calibration matrices
K2 = [7.188560000000e+02 0.000000000000e+00 6.071928000000e+02;...
    0.000000000000e+00 7.188560000000e+02 1.852157000000e+02;...
    0.000000000000e+00 0.000000000000e+00 1.000000000000e+00];
K3 = [ 7.188560000000e+02 0.000000000000e+00 6.071928000000e+02;...
    0.000000000000e+00 7.188560000000e+02 1.852157000000e+02;...
    0.000000000000e+00 0.000000000000e+00 1.000000000000e+00];

% homogeneous transforms between camera 2 and camera3
% camera 2 is left, camera 3 right
H2 = [eye(3) [-0.06 -1.65 0]'];
H3 = [eye(3) [0.48 -1.65 0]'];

%% initialize feature matcher
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

num_images = min(length(images2),length(images3));

%% run VSLAM
debug_flag = true;

for i = 1:1%num_images
    %% get next image sequence
    im_L = imread([img2Path images2(i).name]);
    im_R = imread([img3Path images3(i).name]);

    %% push stereo images into feature detector
    matcherMex('push',rgb2gray(im_L),rgb2gray(im_R));

    %% match features
    matcherMex('match',1);
    matched_features = matcherMex('get_matches',1);
    matched_indices = matcherMex('get_indices',1);
    
    %% extract matched features
    feature_L = matched_features(1:2,:)';
    feature_R = matched_features(3:4,:)';
    
    %% remove outliers using RANSAC
    [x_l, x_r, idx] = GetInliersRANSAC(feature_L,feature_R,0.01,1000);
    
    %% perform linear triangulation of points
    X = LinearTriangulation(K2,K3,H2(:,4),H2(:,1:3),H3(:,4),H3(:,1:3),x_l,x_r);
    if debug_flag && ~isempty(idx)
        %% plot features
        %{d
        figure(1)
        clf
        imshow(im_L);
        hold on
        plot(x_l(:,1),x_l(:,2),'r*')

        figure(2)
        clf
        imshow(im_R);
        hold on
        plot(x_r(:,1),x_r(:,2),'b*')
        %}
        figure(3)
        clf
        mask = sqrt(sum(X.^2,2))<30;
        R_w = [0 0 1; -1 0 0; 0 -1 0];
        showPointCloud(X(mask,:))
        xlabel('x')
        ylabel('y')
        zlabel('z')
    end
    %pause
end

%% close matcher
matcherMex('close');