clear all
close all
addpath('imu')
addpath('cam')
addpath('ref')
addpath('vicon')

ukf_state = [];
vicon_state = [];
acc_meas = [];
vicon_meas = [];

%% declare dataset
dataset = 1;
imu_file = ['imuRaw' num2str(dataset)];
cam_file = ['cam' num2str(dataset)];
vicon_file = ['viconRot' num2str(dataset)];

%% load dataset
load('imu_params.mat')
imu = load(imu_file);
vicon = load(vicon_file);
cam = load(cam_file);
cam.cam(1:floor(end/2),:,:,1) = 0;
%imshow(cam.cam(:,:,:,1))
%% plots
lims = [-1.1 1.1];

p_v = figure(1);
clf
r_vicon = eye(3);
p_vicon.x = plot3([0 r_vicon(1,1)],[0 r_vicon(2,1)],[0 r_vicon(3,1)],'r-*');
hold on
axis equal
grid on
title('Vicon Attitude')
xlabel('x')
ylabel('y')
zlabel('z')
xlim(lims);
ylim(lims);
zlim(lims);
p_vicon.y = plot3([0 r_vicon(1,2)],[0 r_vicon(2,2)],[0 r_vicon(3,2)],'g-*');
p_vicon.z = plot3([0 r_vicon(1,3)],[0 r_vicon(2,3)],[0 r_vicon(3,3)],'b-*');

p_g = figure(2);
clf
r_gyro = eye(3);
p_gyro.x = plot3([0 r_gyro(1,1)],[0 r_gyro(2,1)],[0 r_gyro(3,1)],'r-*');
hold on
axis equal
grid on
title('UKF Attitude')
xlabel('x')
ylabel('y')
zlabel('z')
xlim(lims);
ylim(lims);
zlim(lims);
p_gyro.y = plot3([0 r_gyro(1,2)],[0 r_gyro(2,2)],[0 r_gyro(3,2)],'g-*');
p_gyro.z = plot3([0 r_gyro(1,3)],[0 r_gyro(2,3)],[0 r_gyro(3,3)],'b-*');

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
    1.0097*10^3;...
    1.0099*10^3;...
    1.0099*10^3;...
    10;...
    50;...
    50 ...
    ]);

R = diag(...
    [
    1.5688e1*1.5;...
    4.7045e1*1.5;...
    1.2642e1*1.5;...
    8.4178*10^-5;...
    7.6938*10^-5;...
    1.1999*10^-4 ...
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

for i = 2:length(imu.ts)
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
    
    X = [quatmultiply(q_k',q_W')';...
        bsxfun(@plus,omega_k,W(4:6,:))];
    
    % transform sigma points
    alpha_delta = sqrt(sum(X(5:7,:).^2))*dt;
    e_delta = bsxfun(@rdivide,X(5:7,:)*dt,alpha_delta);
    e_delta(:,alpha_delta == 0) = repmat([0 0 0]',1,sum(alpha_delta==0));
    q_delta = [cos(alpha_delta/2);...
        bsxfun(@times,e_delta,sin(alpha_delta/2))];
    q_k1 = quatmultiply(q_k',q_delta')';
    
    Y = [quatmultiply(X(1:4,:)',q_delta')';...
        X(5:7,:)];
    
    % calculate sigma point mean and covariance
    [~,~,V] = svd((Y(1:4,:)*Y(1:4,:)')/(2*n));
    x_k_mean = [V(:,1)/norm(V(:,1));mean(Y(end-2:end,:),2)];
    
    r_prime = quatmultiply(Y(1:4,:)',x_k_mean(1:4)'.*[1 -1 -1 -1])';
    omega_prime = bsxfun(@minus,Y(5:7,:),x_k_mean(5:7));
    Y_mean_centered = [r_prime; omega_prime];
    W_prime = [bsxfun(@times,2*acos(Y_mean_centered(1,:)) ...
                            ./sin(acos(Y_mean_centered(1,:))),Y_mean_centered(2:4,:));...
            Y_mean_centered(end-2:end,:)];
    P_k_mean = W_prime*W_prime'/(2*n);
    
    % compute a-priori estimate
    Z = [quatrotate(Y(1:4,:)',[0 0 1])';Y(end-2:end,:)];
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
    
    %combine = [x_q; x_k_mean(5:7)]+v_kc;
    
    alpha_v = norm(v_kc(1:3));
    e_v = v_kc(1:3)/norm(v_kc(1:3));
    e_v(:,alpha_v == 0) = repmat([0 0 0]',1,sum(alpha_v==0));
    q_v = [cos(alpha_v/2);...
        e_v*sin(alpha_v/2)];
    q_combined = quatmultiply(x_k_mean(1:4)',q_v')';
    x_k = [q_combined; x_k_mean(5:7)+v_kc(4:6)];
    P_k = P_k_mean-K_k*P_vv*K_k';
    
    
    %% image stitching
    [~,im_idx] = min(abs(imu.ts(i)-cam.ts));
    
    [~,vicon_idx] = min(abs(imu.ts(i)-vicon.ts));
    q_Vicon = rot_to_quat(vicon.rots(:,:,vicon_idx));
    rotated_coord = quatrotate((x_k(1:4).*[1 -1 -1 -1]')',space_coords);
    rotated_coord = quatrotate((q_Vicon.*[1 -1 -1 -1]')',space_coords);
    coord = [asin(rotated_coord(:,3)./space_coord_norm) atan2(rotated_coord(:,2),rotated_coord(:,1))];
    shifted_coord = bsxfun(@plus,coord,[pi/2 pi]);
    idx_coord = ceil(shifted_coord/d_angle);
    pixels = reshape(flipud(imrotate(cam.cam(:,:,:,im_idx),90)), 320*240*3,1);
    lin_coords_r = sub2ind(size(canvas(:,:,1)),idx_coord(:,1),idx_coord(:,2));
    lin_coords_g = lin_coords_r + numel(canvas(:,:,1));
    lin_coords_b = lin_coords_r + numel(canvas(:,:,1))*2;
    canvas([lin_coords_r; lin_coords_g; lin_coords_b]) = pixels;

    %canvas(idx_coord(:,1),idx_coord(:,2),:) = reshape(cam.cam(:,:,:,im_idx), 320*240,3);
    set(im,'cdata',flipud(canvas))
    %{
    %% plotting
    % plot propagated rotation
    rot1 = quat_to_rot(x_k_mean(1:4));
    r_gyro = rot1*eye(3);
    
    set(p_gyro.x,'xdata',[0 r_gyro(1,1)],'ydata',[0 r_gyro(2,1)],'zdata',[0 r_gyro(3,1)]);
    set(p_gyro.y,'xdata',[0 r_gyro(1,2)],'ydata',[0 r_gyro(2,2)],'zdata',[0 r_gyro(3,2)]);
    set(p_gyro.z,'xdata',[0 r_gyro(1,3)],'ydata',[0 r_gyro(2,3)],'zdata',[0 r_gyro(3,3)]);
    set(p_g,'Name',num2str(imu.ts(i)-imu.ts(1)));
    drawnow
    
    [~,vicon_idx] = min(abs(imu.ts(i)-vicon.ts));
    r_vicon = vicon.rots(:,:,vicon_idx)*eye(3);
    set(p_vicon.x,'xdata',[0 r_vicon(1,1)],'ydata',[0 r_vicon(2,1)],'zdata',[0 r_vicon(3,1)]);
    set(p_vicon.y,'xdata',[0 r_vicon(1,2)],'ydata',[0 r_vicon(2,2)],'zdata',[0 r_vicon(3,2)]);
    set(p_vicon.z,'xdata',[0 r_vicon(1,3)],'ydata',[0 r_vicon(2,3)],'zdata',[0 r_vicon(3,3)]);
    set(p_v,'Name',num2str(vicon.ts(vicon_idx)-vicon.ts(1)));
    drawnow

    ukf_state = [ukf_state x_k(1:4)];
    vicon_state = [vicon_state rot_to_quat(vicon.rots(:,:,vicon_idx))];
    
    acc_meas = [acc_meas quatrotate(rot_to_quat(vicon.rots(:,:,vicon_idx))',((imu.vals(1:3,i)-a_b).*params.sf_a.*[-1;-1;1])')'];
    vicon_meas = [vicon_meas quatrotate(rot_to_quat(vicon.rots(:,:,vicon_idx))',[0 0 1])'];
    %}
    stop = toc;
    %fprintf('time: %6.6f\n', imu.ts(i)-imu.ts(1))
    pause(dt-(stop-start))
    start = stop;
end

%{
%%
figure(5)
clf
plot(vicon_state(1,:),'r-')
hold on
plot(ukf_state(1,:),'b-')
grid on
legend('vicon','ukf')

figure(6)
clf
plot(vicon_state(2,:),'r-')
hold on
plot(ukf_state(2,:),'b-')
grid on
legend('vicon','ukf')

figure(7)
clf
plot(vicon_state(3,:),'r-')
hold on
plot(ukf_state(3,:),'b-')
grid on
legend('vicon','ukf')

figure(8)
clf
plot(vicon_state(4,:),'r-')
hold on
plot(ukf_state(4,:),'b-')
grid on
legend('vicon','ukf')
%}