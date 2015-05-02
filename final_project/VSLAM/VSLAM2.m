clear all
close all

%% load data
img2Path = '../dataset/sequences/02/image_2/';
imgType = '*.png';
images2 = dir([img2Path imgType]);

img3Path = '../dataset/sequences/02/image_3/';
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

%{
im2_prev = imread([img2Path images2(1).name]);
oldPoints2 = corner(rgb2gray(im2_prev));
pointTracker2 = vision.PointTracker('MaxBidirectionalError', 2);
initialize(pointTracker2, oldPoints2, rgb2gray(im2_prev));
%}
%{
f = figure(1);
clf
imshow(rgb2gray(im2_prev));
%}
figure(1);
pc = plot3(0,0,0,'k*');%showPointCloud(zeros(3,3));
axis equal
grid on
xlabel('x')
ylabel('y')
zlabel('z')

figure(2);
im2_plot = imshow(rgb2gray(imread([img2Path images2(1).name])));
hold on
im2_feat = plot(0,0,'r+');

figure(3);
im3_plot = imshow(rgb2gray(imread([img3Path images3(1).name])));
hold on
im3_feat = plot(0,0,'b+');
%}
X_full = [];
for i = 2:length(images2)
    
    %% read image sequence
    im2_rgb = imread([img2Path images2(i).name]);
    im3_rgb = imread([img3Path images3(i).name]);
    im2 = rgb2gray(im2_rgb);
    im3 = rgb2gray(im3_rgb);
    
    matcherMex('push',im2);
    matcherMex('match',0);
    p_matched{i} = matcherMex('get_matches',0);
    i_matched{i} = matcherMex('get_indices',0);
    
    %% create features
    points_2 = detectHarrisFeatures(im2);
    points_3 = detectHarrisFeatures(im3);
    [feat_2, valid_points_2] = extractFeatures(im2,points_2);
    [feat_3, valid_points_3] = extractFeatures(im3,points_3);
    
    %% match features
    indexPairs = matchFeatures(feat_2,feat_3);
    feature_2 = double(valid_points_2(indexPairs(:,1),:).Location);
    feature_3 = double(valid_points_3(indexPairs(:,2),:).Location);
    
    set(im2_plot,'cdata',im2)
    set(im2_feat,'xdata',feature_2(:,1),'ydata', feature_2(:,2));
    
    set(im3_plot,'cdata',im3)
    set(im3_feat,'xdata',feature_3(:,1),'ydata', feature_3(:,2));
    %{d
    %% RANSAC outlier rejection
    [x1, x2, idx] = GetInliersRANSAC(feature_2,feature_3,0.5,10000);

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

    %rotated_X = X0*[0 -1 0; 0 0 -1; 1 0 0];

    X_new = NonlinearTriangulation(K3, zeros(3,1), eye(3), C, R, x1, x2, X0);
    rotated_X = X_new*[0 -1 0; 0 0 -1; 1 0 0];
    set(pc,'XData',rotated_X(:,1),'YData',rotated_X(:,2),'ZData',rotated_X(:,3))
    %{
    figure
    showPointCloud(X_new*[0 -1 0; 0 0 -1; 1 0 0])
    axis equal
    xlabel('x')
    ylabel('y')
    zlabel('z')
    fprintf('Time to triangulate: %6.6f seconds\n', t_e-t_s);
    %}
    %{d
    figure(6)
    clf
    imshow(im2)
    hold on
    for j = 1:length(x1)
        plot(x1(j,1), x1(j,2),'b+');
        plot(x2(j,1),x2(j,2),'g+');
        plot([x1(j,1) x2(j,1)],[x1(j,2) x2(j,2)],'r-');
    end
    %}
    drawnow
    %{
    %% assign previous features and camera pose
    [~,vol] = convhull(double(feature_2));
    if vol <= numel(im2)*0.8
        oldPoints2 = corner(im2);
    else
        oldPoints2 = feature_2;
    end
    setPoints(pointTracker2, oldPoints2);
    %}
    
    %C0 = C;
    %R0 = R;
    pause
end