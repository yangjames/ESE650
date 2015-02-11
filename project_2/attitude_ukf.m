clear all
addpath('imu')
addpath('cam')
addpath('ref')
addpath('vicon')

%% declare dataset
dataset = 1;
imu_file = ['imuRaw' num2str(dataset)];
cam_file = ['cam' num2str(dataset)];
vicon_file = ['viconRot' num2str(dataset)];

%% load dataset
load('imu_params.mat')
imu = load(imu_file);
vicon = load(vicon_file);

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

p_a = figure(3);
clf
r_acc = eye(3);
p_acc.x = plot3([0 r_acc(1,1)],[0 r_acc(2,1)],[0 r_acc(3,1)],'r-*');
hold on
axis equal
grid on
title('Accelerometer Attitude')
xlabel('x')
ylabel('y')
zlabel('z')
xlim(lims);
ylim(lims);
zlim(lims);
p_acc.y = plot3([0 r_acc(1,2)],[0 r_acc(2,2)],[0 r_acc(3,2)],'g-*');
p_acc.z = plot3([0 r_acc(1,3)],[0 r_acc(2,3)],[0 r_acc(3,3)],'b-*');

%% constants
% bias, scale factor, high-pass filter
omega_b = [params.roll_bias;params.pitch_bias;params.yaw_bias];
a_b = [params.ax_bias; params.ay_bias; params.az_bias];
a_tol = 0.7;
cutoff_high = 0.001;
RC_high = 1/(cutoff_high*2*pi);

% process noise
Q = eye(6);

%% initialize
% initial state mean and covariance
q_k = [1 0 0 0]'; 
x_k = [q_k; zeros(3,1)];
P_k = eye(6);
n = length(P_k);
omega_k = zeros(3,1);

Y = zeros(7,12);
Z = zeros(7,12);
H_2 = zeros(4,12);
H_1 = zeros(3,12);

vicon_idx = 1;
tic
start = toc;
prev = zeros(7,1);
for i = 2:length(imu.ts)
    % get dt and high-passed angular velocity vector
    dt = imu.ts(i)-imu.ts(i-1);
    alpha_hp = RC_high/(dt+RC_high);
    omega_1 = (imu.vals([5 6 4],i)-omega_b)*params.sf_w;
    omega_0 = (imu.vals([5 6 4],i-1)-omega_b)*params.sf_w;
    omega_k = omega_k*alpha_hp + alpha_hp*(omega_1-omega_0);
    
    % obtain sigma points using Cholesky decomposition
    S = chol(P_k+Q,'upper');
    W = [-sqrt(2*n)*S sqrt(2*n)*S];
    
    % obtain sigma points using diagonalization
    %[~,S,V] = svd(P_k+Q);
    %W = [sqrt(2*n)*V*sqrt(S)*V' -sqrt(2*n)*V*sqrt(S)*V']
    
    % update sigma points
    alpha_W = sqrt(sum(W(1:3,:).^2));
    
    e_W = bsxfun(@rdivide,W(1:3,:),alpha_W);
    e_W(:,alpha_W == 0) = repmat([1 0 0]',1,sum(alpha_W==0));
    q_W = [cos(alpha_W/2);...
        bsxfun(@times,e_W,sin(alpha_W/2))];
    X = [bsxfun(@quat_mult,q_k,q_W);...
        bsxfun(@plus,omega_k,W(4:6,:))];
    
    % propagate state
    alpha_delta = norm(omega_k)*dt;
    e_delta = omega_k/norm(omega_k);
    q_delta = [cos(alpha_delta/2);e_delta*sin(alpha_delta/2)];
    q_k1 = quat_mult(q_k,q_delta);
    Y = [bsxfun(@quat_mult,X(1:4,:),q_delta);...
        X(5:7,:)];
    
    % calculate sigma point mean and covariance
    [~,~,V] = svd(Y(1:4,:)*Y(1:4,:)'/size(Y,2));
    x_k_mean = [V(:,1);mean(Y(end-2:end,:),2)];
    
    Y_mean_centered = bsxfun(@minus,Y,x_k_mean);
    W_prime = [bsxfun(@times,2*acos(Y_mean_centered(1,:))./sin(acos(Y_mean_centered(1,:))),Y(2:4,:));...
        Y_mean_centered(end-2:end,:)];
    P_k_mean = W_prime*W_prime'/size(W_prime,2);
    
    %{
    % update with accelerometer data if magnitude is reasonable
    a = (imu.vals(1:3,i)-a_b)*params.sf_a.*[-1;-1;1];
    if norm(a) < 9.81+a_tol && norm(a) > 9.81-a_tol ...
            && abs(a(1)) < 9.81-a_tol ...
            && abs(a(2)) < 9.81-a_tol
        temp = bsxfun(@quat_mult,Y(1:4,:),[0;a]);
        for j = 1:size(Y,2)
            H_2(:,i) = quat_mult(temp(:,i),Y(1:4,i).*[1;-1;-1;-1]);
        end
    end
    %}
    
    q_k = x_k_mean(1:4);
    %% plotting
    % plot propagated rotation
    rot1 = quat_to_rot(x_k_mean(1:4));
    r_gyro = rot1*eye(3);
    
    set(p_gyro.x,'xdata',[0 r_gyro(1,1)],'ydata',[0 r_gyro(2,1)],'zdata',[0 r_gyro(3,1)]);
    set(p_gyro.y,'xdata',[0 r_gyro(1,2)],'ydata',[0 r_gyro(2,2)],'zdata',[0 r_gyro(3,2)]);
    set(p_gyro.z,'xdata',[0 r_gyro(1,3)],'ydata',[0 r_gyro(2,3)],'zdata',[0 r_gyro(3,3)]);
    set(p_g,'Name',num2str(imu.ts(i)-imu.ts(1)));
    drawnow
    
    lt = vicon.ts-vicon.ts(1) < imu.ts(i)-imu.ts(1);
    t_stamp = vicon.ts(lt);
    r_vicon = vicon.rots(:,:,length(t_stamp))*eye(3);
    set(p_vicon.x,'xdata',[0 r_vicon(1,1)],'ydata',[0 r_vicon(2,1)],'zdata',[0 r_vicon(3,1)]);
    set(p_vicon.y,'xdata',[0 r_vicon(1,2)],'ydata',[0 r_vicon(2,2)],'zdata',[0 r_vicon(3,2)]);
    set(p_vicon.z,'xdata',[0 r_vicon(1,3)],'ydata',[0 r_vicon(2,3)],'zdata',[0 r_vicon(3,3)]);
    set(p_v,'Name',num2str(vicon.ts(length(t_stamp))-vicon.ts(1)));
    drawnow

    stop = toc;
    %fprintf('time: %6.6f\n', imu.ts(i)-imu.ts(1))
    pause(dt-(stop-start))
    start = stop;

end