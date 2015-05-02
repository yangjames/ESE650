clear all
close all
addpath('../libviso2/matlab/')
%% load data
img2Path = '../dataset/sequences/05/image_2/';
imgType = '*.png';
images2 = dir([img2Path imgType]);

img3Path = '../dataset/sequences/05/image_3/';
images3 = dir([img3Path imgType]);

K2 = [7.188560000000e+02 0.000000000000e+00 6.071928000000e+02;...
    0.000000000000e+00 7.188560000000e+02 1.852157000000e+02;...
    0.000000000000e+00 0.000000000000e+00 1.000000000000e+00];
K3 = [ 7.188560000000e+02 0.000000000000e+00 6.071928000000e+02;...
    0.000000000000e+00 7.188560000000e+02 1.852157000000e+02;...
    0.000000000000e+00 0.000000000000e+00 1.000000000000e+00];

H2 = [eye(3) [0 0 1.65]'];
H3 = [eye(3) [0 0.54 1.65]'];

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
R0 = eye(3);
C0 = zeros(3,1);

% initialize feature tracker
im2_prev = imread([img2Path images2(1).name]);
matcherMex('push',rgb2gray(im2_prev));


X_full = [];
for i = 2:length(images2)
    %% read image sequence
    im2_rgb = imread([img2Path images2(i).name]);
    im2 = rgb2gray(im2_rgb);
    
    %% create and match features
    matcherMex('push',im2);
    matcherMex('match',0);
    features = matcherMex('get_matches',0);
    matched_indices = matcherMex('get_indices',0);
    
    feature_1 = features(1:2,:)';
    feature_2 = features(3:4,:)';
    %{d
    %% RANSAC outlier rejection
    [x1, x2, idx] = GetInliersRANSAC(feature_1,feature_2,1,1000);
    %x1 = feature_1;
    %x2 = feature_2;

    %% fundamental and essential matrix estimation
    F = EstimateFundamentalMatrix(x1,x2);
    E = EssentialMatrixFromFundamentalMatrix(F,K3);
    
    %% camera pose extraction and feature triangulation
    [Cset, Rset] = ExtractCameraPose(E);
    Xset = cell(4,1);
    for j = 1:4
        Xset{j} = LinearTriangulation(K3, C0, R0, Cset{j}, Rset{j}, x1, x2);
    end
    
    %% camera pose disambiguation
    [C, R, X0] = DisambiguateCameraPose(Cset,Rset,Xset);
    %X_new = X0;
    X_new = NonlinearTriangulation(K3, zeros(3,1), eye(3), C, R, x1, x2, X0);
    rotated_X = X_new*[0 -1 0; 0 0 -1; 1 0 0];
    mask_x = sqrt(sum(rotated_X.^2,2))<30;
    set(pc,'XData',rotated_X(mask_x,1),'YData',rotated_X(mask_x,2),'ZData',rotated_X(mask_x,3))
    %{
    figure
    showPointCloud(X_new*[0 -1 0; 0 0 -1; 1 0 0])
    axis equal
    xlabel('x')
    ylabel('y')
    zlabel('z')
    %}
    %{d
    figure(2)
    clf
    imshow(im2)
    hold on
    for j = 1:length(x1)
        plot(x1(j,1), x1(j,2),'b+');
        plot(x2(j,1), x2(j,2),'g+');
        plot([x1(j,1) x2(j,1)],[x1(j,2) x2(j,2)],'r-');
    end
    %}
    drawnow
    
    %C0 = C;
    %R0 = R;
end