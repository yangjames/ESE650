clear all
close all
addpath('imu')
addpath('cam')
addpath('ref')
addpath('vicon')
addpath('P2_TEST')

ukf_state = [];
vicon_state = [];
acc_meas = [];
vicon_meas = [];

%% declare dataset
dataset = 10;
imu_file = ['imu_test'];
cam_file = ['cam_test'];
%imu_file= ['imuRaw' num2str(dataset)];
%cam_file = ['cam' num2str(dataset)];

%% load dataset
load('imu_params.mat')
imu = load(imu_file);
cam = load(cam_file);

%% constants
% bias, scale factor
omega_b = [params.roll_bias;params.pitch_bias;params.yaw_bias];
a_b = [params.ax_bias; params.ay_bias; params.az_bias];

%% initialize
% initial state mean and covariance
x_k = [[1 0 0 0]'; zeros(3,1)];
P_k = 1e-2*eye(6);
n = length(P_k);

% process and sensor noise
Q = diag(...
    [
    78.98;...
    78.98;...
    78.98;...
    18.9;...
    18.9;...
    18.9 ...
    ]);

R = diag(...
    [
    0.5;...
    0.5;...
    0.5;...
    0.001;...
    0.001;...
    0.001 ...
    ]);

% plotting stuff
vicon_idx = 1;
tic
start = toc;
prev = zeros(7,1);

num_pix = 500;
d_angle = pi/num_pix;

canvas = uint8(zeros(num_pix,2*num_pix,3));
figure(4)
im = imshow(canvas);
[row,col] = ind2sub([320,240],1:320*240);
pixel_coords = [row; col];
space_coords = pixel_to_world(pixel_coords)';
space_coord_norm = sqrt(sum(space_coords.^2,2));
prev_idx = 0;
Y = zeros(7,13);
X = Y;
for i = 2:length(imu.ts)-700
    % get dt
    dt = imu.ts(i)-imu.ts(i-1);
    q_k = x_k(1:4);
    omega_k = x_k(end-2:end);
    
    % obtain sigma points using Cholesky decomposition
    S = chol(P_k+Q);
    W = [sqrt(n)*S zeros(6,1) -sqrt(n)*S];
    
    % propagate sigma points
    alpha_W = sqrt(sum(W(1:3,:).^2));
    e_W = bsxfun(@rdivide,W(1:3,:),alpha_W);
    e_W(:,alpha_W == 0) = repmat([0 0 0]',1,sum(alpha_W==0));
    q_W = [cos(alpha_W/2);...
        bsxfun(@times,e_W,sin(alpha_W/2))];
    
    for j = 1:size(q_W,2)
        X(:,j) = [quat_mult(q_k,q_W(:,j)); omega_k+W(4:6,j)];
    end
    X(1:4,:) = bsxfun(@rdivide,X(1:4,:),sqrt(sum(X(1:4,:).^2)));
    
    % transform sigma points
    alpha_delta = sqrt(sum(X(5:7,:).^2))*dt;
    e_delta = bsxfun(@rdivide,X(5:7,:)*dt,alpha_delta);
    e_delta(:,alpha_delta == 0) = repmat([0 0 0]',1,sum(alpha_delta==0));
    q_delta = [cos(alpha_delta/2);...
        bsxfun(@times,e_delta,sin(alpha_delta/2))];
    q_k1 = bsxfun(@quat_mult,q_k,q_delta);
    
    for j = 1:size(q_delta,2)
        Y(:,j) = [quat_mult(X(1:4,j),q_delta(:,j)); X(5:7,j)];
    end
    
    Y(1:4,:) = bsxfun(@rdivide,Y(1:4,:),sqrt(sum(Y(1:4,:).^2)));
    
    % calculate sigma point mean and covariance
    [~,~,V] = svd((Y(1:4,:)*Y(1:4,:)')/(2*n));
    x_k_mean = [V(:,1)/norm(V(:,1));mean(Y(end-2:end,:),2)];
    
    r_prime = bsxfun(@quat_mult,quat_conj(x_k_mean(1:4)),Y(1:4,:));
    omega_prime = bsxfun(@minus,Y(5:7,:),x_k_mean(5:7));
    Y_mean_centered = [r_prime; omega_prime];
    W_prime = [bsxfun(@times,2*acos(Y_mean_centered(1,:)) ...
                            ./sin(acos(Y_mean_centered(1,:))),Y_mean_centered(2:4,:));...
            Y_mean_centered(end-2:end,:)];
    P_k_mean = W_prime*W_prime'/(2*n);
    
    % compute a-priori estimate
    g_rot = bsxfun(@quat_mult,quat_conj(Y(1:4,:)),[0 0 0 1]');
    for j = 1:size(Y,2)
        g_rot(:,j) = quat_mult(g_rot(:,j),Y(1:4,j));
    end
    
    Z = [g_rot(2:4,:);Y(end-2:end,:)];
    z_k_mean = mean(Z,2);
    
    Z_mean_centered = bsxfun(@minus,Z,z_k_mean);
    P_zz = (Z_mean_centered*Z_mean_centered')/(2*n);
    P_vv = P_zz + R;
    
    P_xz = (W_prime*Z_mean_centered')/(2*n);
    
    % calculate Kalman gain, residual, and state update
    K_k = P_xz/P_vv;
    v_k = [(imu.vals(1:3,i)-a_b).*params.sf_a.*[-1;-1;1];...
        (imu.vals([5 6 4],i)-omega_b).*params.sf_w]-z_k_mean;
    v_kc = K_k*v_k;
    
    x_q = quat_to_vec(x_k_mean(1:4)');
    
    alpha_v = norm(v_kc(1:3));
    e_v = v_kc(1:3)/norm(v_kc(1:3));
    e_v(:,alpha_v == 0) = repmat([0 0 0]',1,sum(alpha_v==0));
    q_v = [cos(alpha_v/2);...
        e_v*sin(alpha_v/2)];
    q_combined = bsxfun(@quat_mult, x_k_mean(1:4),q_v);
    x_k = [bsxfun(@rdivide,q_combined,sqrt(sum(q_combined.^2))); x_k_mean(5:7)+v_kc(4:6)];
    P_k = P_k_mean-K_k*P_vv*K_k';
    
    
    %% image stitching
    [~,im_idx] = min(abs(imu.ts(i)-cam.ts));
    if prev_idx ~= im_idx
        rotated_coord = quatrotate((x_k(1:4).*[1 -1 -1 -1]')',space_coords);
        coord = [asin(rotated_coord(:,3)./space_coord_norm) atan2(rotated_coord(:,2),rotated_coord(:,1))];
        shifted_coord = bsxfun(@plus,coord,[pi/2 pi]);
        idx_coord = ceil(shifted_coord/d_angle);
        pixels = reshape(flipud(imrotate(cam.cam(:,:,:,im_idx),90)), 320*240*3,1);
        lin_coords_r = sub2ind(size(canvas(:,:,1)),idx_coord(:,1),idx_coord(:,2));
        lin_coords_g = lin_coords_r + numel(canvas(:,:,1));
        lin_coords_b = lin_coords_r + numel(canvas(:,:,1))*2;
        canvas([lin_coords_r; lin_coords_g; lin_coords_b]) = pixels;
        prev_idx = im_idx;
    end

    set(im,'cdata',rot90(canvas,2))
    stop = toc;
    pause(dt-(stop-start))
    start = stop;
end