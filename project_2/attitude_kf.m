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
figure(1)
clf
r_vicon = eye(3);
p_vicon.x = plot3([0 r_vicon(1,1)],[0 r_vicon(2,1)],[0 r_vicon(3,1)],'r-*');
hold on
axis equal
grid on
xlim([-1 1]);
ylim([-1 1]);
zlim([-1 1]);
p_vicon.y = plot3([0 r_vicon(1,2)],[0 r_vicon(2,2)],[0 r_vicon(3,2)],'g-*');
p_vicon.z = plot3([0 r_vicon(1,3)],[0 r_vicon(2,3)],[0 r_vicon(3,3)],'b-*');

figure(2)
clf
r_gyro = eye(3);
p_gyro.x = plot3([0 r_gyro(1,1)],[0 r_gyro(2,1)],[0 r_gyro(3,1)],'r-*');
hold on
axis equal
grid on
xlim([-1 1]);
ylim([-1 1]);
zlim([-1 1]);
p_gyro.y = plot3([0 r_gyro(1,2)],[0 r_gyro(2,2)],[0 r_gyro(3,2)],'g-*');
p_gyro.z = plot3([0 r_gyro(1,3)],[0 r_gyro(2,3)],[0 r_gyro(3,3)],'b-*');

figure(3)
clf
r_acc = eye(3);
p_acc.x = plot3([0 r_acc(1,1)],[0 r_acc(2,1)],[0 r_acc(3,1)],'r-*');
hold on
axis equal
grid on
xlim([-1 1]);
ylim([-1 1]);
zlim([-1 1]);
p_acc.y = plot3([0 r_acc(1,2)],[0 r_acc(2,2)],[0 r_acc(3,2)],'g-*');
p_acc.z = plot3([0 r_acc(1,3)],[0 r_acc(2,3)],[0 r_acc(3,3)],'b-*');

%% initialize
Q = zeros(6,1);%[0 0 0 0 0.1 0.1 0.1]'; % dynamics noise
R = [0.1 0.1 0.1]'; % accelerometer noise
w_b = [params.roll_bias;params.pitch_bias;params.yaw_bias];
a_b = [params.ax_bias; params.ay_bias; params.az_bias];
a_tol = 1;
cutoff_high = 0.001;
RC_high = 1/(cutoff_high*2*pi);

att_0 = [0 0 0]';
w_0 = (imu.vals([5 6 4],1) - w_b)*params.sf_w;
mu_0 = [att_0; w_0];
S_0 = eye(length(mu_0));
vicon_idx = 1;
tic
start = toc;
w=[0;0;0];
for i = 2:length(imu.ts)
     % get dt and high-passed angular velocity vector
    dt = imu.ts(i)-imu.ts(i-1);
    alpha_hp = RC_high/(dt+RC_high);
    w_1 = (imu.vals([5 6 4],i)-w_b)*params.sf_w;
    w_0 = (imu.vals([5 6 4],i-1)-w_b)*params.sf_w;
    w = w*alpha_hp + alpha_hp*(w_1-w_0);

    % propagate first order integration
    A = [eye(3) dt*eye(3);...
        zeros(3) eye(3)];
    mu = A*[mu_0(1:3);w];
    S = A*S_0*A'+diag(Q);
    
    % update with accelerometer data if magnitude is reasonable
    a = (imu.vals(1:3,i)-a_b)*params.sf_a.*[-1;-1;1];
    if norm(a) < 9.81+a_tol && norm(a) > 9.81-a_tol
        roll = (atan2(a(2), a(3)))*180/pi;
        pitch = (atan2(-a(1), a(3)))*180/pi;
        if roll < -180
            roll = roll + 180;
        elseif roll >= 180
            roll = roll - 180;
        end
        if pitch < -180
            pitch = pitch + 180;
        elseif pitch >= 180
            pitch = pitch - 180;
        end
    end
    
    % assign mu and sigma for next time step
    mu_0 = mu;
    S_0 = S;
    
    %% plotting
    % get rotation matrices
    rot1 = euler_to_rot(mu_0(1),mu_0(2),mu_0(3));
    rot2 = euler_to_rot(roll*pi/180,pitch*pi/180,mu(3));
    
    % plot rotations
    r_gyro = rot1*eye(3);
    set(p_gyro.x,'xdata',[0 r_gyro(1,1)],'ydata',[0 r_gyro(2,1)],'zdata',[0 r_gyro(3,1)]);
    set(p_gyro.y,'xdata',[0 r_gyro(1,2)],'ydata',[0 r_gyro(2,2)],'zdata',[0 r_gyro(3,2)]);
    set(p_gyro.z,'xdata',[0 r_gyro(1,3)],'ydata',[0 r_gyro(2,3)],'zdata',[0 r_gyro(3,3)]);
    drawnow
    
    r_acc = rot2*eye(3);
    set(p_acc.x,'xdata',[0 r_acc(1,1)],'ydata',[0 r_acc(2,1)],'zdata',[0 r_acc(3,1)]);
    set(p_acc.y,'xdata',[0 r_acc(1,2)],'ydata',[0 r_acc(2,2)],'zdata',[0 r_acc(3,2)]);
    set(p_acc.z,'xdata',[0 r_acc(1,3)],'ydata',[0 r_acc(2,3)],'zdata',[0 r_acc(3,3)]);
    drawnow
    
    if vicon.ts(vicon_idx)-vicon.ts(1) < imu.ts(i)-imu.ts(1)
        if vicon_idx+1 <= length(vicon.ts)
            vicon_idx = vicon_idx+1;
            r_vicon = vicon.rots(:,:,vicon_idx)*eye(3);
            set(p_vicon.x,'xdata',[0 r_vicon(1,1)],'ydata',[0 r_vicon(2,1)],'zdata',[0 r_vicon(3,1)]);
            set(p_vicon.y,'xdata',[0 r_vicon(1,2)],'ydata',[0 r_vicon(2,2)],'zdata',[0 r_vicon(3,2)]);
            set(p_vicon.z,'xdata',[0 r_vicon(1,3)],'ydata',[0 r_vicon(2,3)],'zdata',[0 r_vicon(3,3)]);
            drawnow
        end
    end

    stop = toc;
    pause(dt-(stop-start))
    fprintf('time: %6.6f\n', imu.ts(i)-imu.ts(1))
    start = stop;
end