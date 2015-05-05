clear all
close all

addpath('../libviso2/matlab/')
%vid = VideoWriter('reconstruct_1_1');
%open(vid)
%% load data
imgType = '*.png';
sequence  = '../dataset/sequences/01/';

img2Path = [sequence 'image_2/'];
images2 = dir([img2Path imgType]);

img3Path = [sequence 'image_3/'];
images3 = dir([img3Path imgType]);

% camera 2 and 3 calibration matrices
K = [7.188560000000e+02 0.000000000000e+00 6.071928000000e+02;...
    0.000000000000e+00 7.188560000000e+02 1.852157000000e+02;...
    0.000000000000e+00 0.000000000000e+00 1.000000000000e+00];


%% initialize feature matcher
param.nms_n                  = 2;   % non-max-suppression: min. distance between maxima (in pixels)
param.nms_tau                = 50;  % non-max-suppression: interest point peakiness threshold
param.match_binsize          = 100;  % matching bin width/height (affects efficiency only)
param.match_radius           = 300; % matching radius (du/dv in pixels)
param.match_disp_tolerance   = 1;   % du tolerance for stereo matches (in pixels)
param.outlier_disp_tolerance = 5;   % outlier removal: disparity tolerance (in pixels)
param.outlier_flow_tolerance = 5;   % outlier removal: flow tolerance (in pixels)
param.multi_stage            = 1;   % 0=disabled,1=multistage matching (denser and faster)
param.half_resolution        = 1;   % 0=disabled,1=match at half resolution, refine at full resolution
param.refinement             = 2;   % refinement (0=none,1=pixel,2=subpixel)

% init matcher
matcherMex('close');
matcherMex('init',param);

num_images = min(length(images2),length(images3));

%% run VSLAM
debug_flag = true;

start_point = 1;
% get initial features
im_L_p = imread([img2Path images2(start_point).name]);
im_R_p = imread([img3Path images3(start_point).name]);
matcherMex('push',rgb2gray(im_L_p),rgb2gray(im_R_p));

im_L   = imread([img2Path images2(start_point+1).name]);
im_R   = imread([img3Path images3(start_point+1).name]);
matcherMex('push',rgb2gray(im_L),rgb2gray(im_R));

im_size = size(im_L);
num_pixels = uint64(prod(im_size(1:2)));

% match features
matcherMex('match',2);
matched_features = matcherMex('get_matches',2);

% extract matched features
Mx = matched_features([1 3 5 7],:)';
My = matched_features([2 4 6 8],:)';
n_features = size(matched_features,2);
%{
figure(1)
clf
imshow(im_L_p)
hold on
plot(matched_features(1,:),matched_features(2,:),'g*')
%}
% remove outliers between all images using RANSAC
[x_l_p, x_r_p, idx_p] = GetInliersRANSAC([Mx(:,1) My(:,1)],[Mx(:,2) My(:,2)],0.05,1000);
ReconX(idx_p) = 1;
%{
figure(2)
clf
imshow(im_L_p)
hold on
plot(x_l_p(:,1),x_l_p(:,2),'g*')
drawnow
while(1)
end
%}
% perform linear triangulation of points
R0 = [1 0 0; 0 cos(0.08) -sin(0.08); 0 sin(0.08) cos(0.08)];
R_p_L = R0;
R_p_R = R0;
C_p_L = zeros(3,1);
C_p_R = -[0.54 0 0]';
X_p = LinearTriangulation(K, C_p_L, R_p_L, C_p_R, R_p_R, x_l_p, x_r_p);
%{
figure(1)
clf
imshow(im_L_p)
hold on
P = K*R0*[eye(3) -C_p_L];
x_projected = bsxfun(@rdivide,P(1:2,:)*[X_p ones(sum(idx_p),1)]',P(3,:)*[X_p ones(sum(idx_p),1)]')';
plot(x_projected(:,1),x_projected(:,2),'g*')
plot(x_l_p(:,1),x_l_p(:,2),'*','color',[1 .5 0])
drawnow
while(1)
end
%}
% perfrom nonlinear triangulation of points
X_p = NonlinearTriangulation( K, C_p_L, R_p_L, C_p_R, R_p_R, x_l_p, x_r_p, X_p);
%{
figure(1)
clf
imshow(im_L_p)
hold on
P = K*R0*[eye(3) -C_p_L];
x_projected = bsxfun(@rdivide,P(1:2,:)*[X_p ones(sum(idx_p),1)]',P(3,:)*[X_p ones(sum(idx_p),1)]')';
plot(x_projected(:,1),x_projected(:,2),'g*')
plot(x_l_p(:,1),x_l_p(:,2),'*','color',[1 .5 0])
drawnow
while(1)
end
%}
R_w = [0 1 0; 0 0 1; -1 0 0];
pos = [0 0.06 1.65 0 -0.48 1.65]';
X = [];
C = [];

%% start path and feature plot
figure(1)
%figure('units','normalized','outerposition',[0 0 1 1])
clf
%subplot(3,1,1)
path_plot_L = plot3(0,0,0,'r-');
hold on
path_plot_R = plot3(0,0,0,'b-');
axis equal
grid on
xlabel('x')
ylabel('y')
zlabel('z')

figure(2)
subplot(3,1,2)
im_plot = imshow(im_L);
hold on
feature_plot = plot(0,0,'r*');

for i = start_point+2:num_images
    %% get pose of left and right cameras at next time frame
    x_L = [Mx(idx_p,3) My(idx_p,3)];%feature_L(idx_p,:);
    C_pnp_L=[];
    R_pnp_L=[];
    n_bytes = 0;
    while isempty(C_pnp_L) || isempty(R_pnp_L)
        fprintf(repmat('\b',1,n_bytes))
        n_bytes = fprintf('no pose estimate yet...');
        [C_pnp_L, R_pnp_L] = PnPRANSAC(X_p, x_L, K,0.1,1000);
    end
    disp('Nonlinear PnP');
    [C_L, R_L] = NonlinearPnP(X_p, x_L, K, C_pnp_L, R_pnp_L);
    
    x_R = [Mx(idx_p,4) My(idx_p,4)];%feature_L(idx_p,:);
    C_pnp_R=[];
    R_pnp_R=[];
    n_bytes = 0;
    while isempty(C_pnp_R) || isempty(R_pnp_R)
        fprintf(repmat('\b',1,n_bytes))
        n_bytes = fprintf('no pose estimate yet...');
        [C_pnp_R, R_pnp_R] = PnPRANSAC(X_p, x_R, K,0.1,1000);
    end
    disp('Nonlinear PnP');
    [C_R, R_R] = NonlinearPnP(X_p, x_R, K, C_pnp_R, R_pnp_R);
    %{
    figure(1)
    clf
    imshow(im_L)
    hold on
    P = K*R_pnp_L*[eye(3) -C_pnp_L];
    x_projected = bsxfun(@rdivide,P(1:2,:)*[X_p ones(sum(idx_p),1)]',P(3,:)*[X_p ones(sum(idx_p),1)]')';
    plot(x_projected(:,1),x_projected(:,2),'g*')
    plot(x_L(:,1),x_L(:,2),'*','color',[1 .5 0])
    drawnow
    while(1)
    end
    %}
    %% triangulate points between new frames
    % remove outliers between all images using RANSAC
    [x_l, x_r, idx] = GetInliersRANSAC([Mx(:,3) My(:,3)],[Mx(:,4) My(:,4)],0.05,1000);

    set(feature_plot,'XData',x_l(:,1),'YData',x_l(:,2));
    % perform linear triangulation of points
    X_LR = LinearTriangulation(K, C_L, R_L, C_R, R_R, x_l, x_r);

    % perfrom nonlinear triangulation of points
    X_LR = NonlinearTriangulation(K, C_L, R_L, C_R, R_R, x_l, x_r, X_LR);

    %% bundle adjustment
    X3D = zeros(n_features,3);
    X3D(idx_p,:) = X_p;
    X3D(idx,:) = X_LR;
    
    C3D = uint8(zeros(n_features,3));
    r_linidx = uint64(sub2ind(im_size(1:2), x_l_p(:,2), x_l_p(:,1)));
    g_linidx = r_linidx + num_pixels;
    b_linidx = g_linidx + num_pixels;
    r_val = im_L_p(r_linidx);
    g_val = im_L_p(g_linidx);
    b_val = im_L_p(b_linidx);
    C3D(idx_p,:) = [r_val(:) g_val(:) b_val(:)];
    
    r_linidx = uint64(sub2ind(im_size(1:2), x_l(:,2), x_l(:,1)));
    g_linidx = r_linidx + num_pixels;
    b_linidx = g_linidx + num_pixels;
    r_val = im_L(r_linidx);
    g_val = im_L(g_linidx);
    b_val = im_L(b_linidx);
    C3D(idx,:) = [r_val(:) g_val(:) b_val(:)];
    
    ReconX = idx | idx_p;
    Mx_bundle = matched_features([1 3 5 7],:)';
    My_bundle = matched_features([2 4 6 8],:)';
    Cr_set = {C_p_L, C_p_R, C_L, C_R};
    Rr_set = {R_p_L, R_p_R, R_L, R_R};
    %disp('Bundle adjustment');
    %[Cr_set, Rr_set, X3D] = BundleAdjustment(K, Cr_set, Rr_set, X3D, ReconX, Mx_bundle, My_bundle);

    C_L = Cr_set{3};
    C_R = Cr_set{4};
    R_L = Rr_set{3};
    R_R = Rr_set{4};
    
    %% store position and reconstruction data and plot
    pos = [pos [R_w'*(R_L'*C_L)+[0 0.06 1.65]';R_w'*(R_R*C_R)+[0 0.06 1.65]']];
    set(path_plot_L,'xdata',pos(1,:),'ydata',pos(2,:),'zdata',pos(3,:));
    set(path_plot_R,'xdata',pos(4,:),'ydata',pos(5,:),'zdata',pos(6,:));
    
    X = [X; bsxfun(@plus,X3D(ReconX,:)*R_w,[0 0.06 2*1.65])];
    C = [C; C3D(ReconX,:)];
    figure(3)
    %clf
    %subplot(3,1,3)
    mask = sqrt(sum(bsxfun(@minus,X,pos(1:3,end)').^2,2)) < 30;% & (X(:,3)-1.65 > 0.1);
    showPointCloud(X(mask,:), C(mask,:));
    drawnow
    frame = getframe(gca);
    writeVideo(vid,frame);

    %% get new images and process them
    % match features
    im_L_p = im_L;
    im_R_p = im_R;
    im_L   = imread([img2Path images2(i).name]);
    im_R   = imread([img3Path images3(i).name]);
       
    %{
    matcherMex('close');
    matcherMex('init',param);
    matcherMex('push',rgb2gray(im_L_p), rgb2gray(im_R_p));
    %}
    matcherMex('push',rgb2gray(im_L),rgb2gray(im_R));
    set(im_plot,'CData',im_L);
    

    matcherMex('match',2);
    matched_features = matcherMex('get_matches',2);
    % extract matched features
    Mx = matched_features([1 3 5 7],:)';
    My = matched_features([2 4 6 8],:)';
    n_features = size(matched_features,2);

    % remove outliers between all images using RANSAC
    [x_l_p, x_r_p, idx_p] = GetInliersRANSAC([Mx(:,1) My(:,1)],[Mx(:,2) My(:,2)],0.05,1000);
    ReconX(idx_p) = 1;

    % perform linear triangulation of points
    R_p_L = Rr_set{3};
    R_p_R = Rr_set{4};
    C_p_L = Cr_set{3};
    C_p_R = Cr_set{4};
    X_p = LinearTriangulation(K, C_p_L, R_p_L, C_p_R, R_p_R, x_l_p, x_r_p);

    % perfrom nonlinear triangulation of points
    X_p = NonlinearTriangulation( K, C_p_L, R_p_L, C_p_R, R_p_R, x_l_p, x_r_p, X_p);
end

%% close matcher
fprintf('\n')
close(vid)
matcherMex('close');