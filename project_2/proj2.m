clear all
addpath('imu')
addpath('cam')
addpath('ref')
addpath('vicon')

dataset = 4;
imu_file = ['imuRaw' num2str(dataset)];
cam_file = ['cam' num2str(dataset)];
vicon_file = ['viconRot' num2str(dataset)];

%% load vicon data and plot
vicon = load(vicon_file);
plotvicon = false;
if plotvicon
    figure(1)
    rotplot(vicon.rots(:,:,1));
    xlabel('x')
    ylabel('y')
    zlabel('z')
    for i = 2:length(vicon.ts)
        rotplot(vicon.rots(:,:,i));
        title(num2str(vicon.ts(i)-vicon.ts(1)))
        drawnow
    end
end

%% load imu data and plot
% acceleration data
imu = load(imu_file);
ax_bias = 510;
ay_bias = 501;
az_bias = 502;
sf_a = [0.0097 0.0097 0.0098]';

figure(2)
clf
plot(imu.ts-imu.ts(1),-(imu.vals(1,:)-ax_bias)*sf_a(1),'r-')
hold on
plot(imu.ts-imu.ts(1),-(imu.vals(2,:)-ay_bias)*sf_a(2),'g-')
plot(imu.ts-imu.ts(1),(imu.vals(3,:)-az_bias)*sf_a(3),'b-')
grid on
title('IMU acceleration plot')
legend('x','y','z')
ylabel('acceleration [g]')
xlabel('time [s]')

var(-(imu.vals(1,1:200)-ax_bias)*sf_a(1))
var(-(imu.vals(2,1:200)-ay_bias)*sf_a(2))
var((imu.vals(3,1:200)-az_bias)*sf_a(3))

% gyro data
size(imu.ts-imu.ts(1))
size(imu.vals(5,:))
roll_bias = 373.75;
pitch_bias = 375.6;
yaw_bias = 370.047;
sf_w = [0.018 0.018 0.018]';
figure(3)
clf
plot(imu.ts-imu.ts(1),(imu.vals(5,:)-roll_bias)*sf_w(1),'r-')
hold on
plot(imu.ts-imu.ts(1),(imu.vals(6,:)-pitch_bias)*sf_w(2),'g-')
plot(imu.ts-imu.ts(1),(imu.vals(4,:)-yaw_bias)*sf_w(3),'b-')
grid on
title('IMU gyro plot')
legend('roll','pitch','yaw')
ylabel('angular rate [rad/s]')
xlabel('time [s]')

var((imu.vals(5,1:200)-roll_bias)*sf_w(1))
var((imu.vals(6,1:200)-pitch_bias)*sf_w(2))
var((imu.vals(4,1:200)-yaw_bias)*sf_w(3))

% angular rate integration plot
dt = imu.ts(2:end)-imu.ts(1:end-1);
roll = cumsum((imu.vals(5,2:end)-roll_bias)*sf_w(1).*dt);
pitch = cumsum((imu.vals(6,2:end)-pitch_bias)*sf_w(2).*dt);
yaw = cumsum((imu.vals(4,2:end)-yaw_bias)*sf_w(3).*dt);
figure(4)
clf
dt = imu.ts(2:end)-imu.ts(1);

plot(dt,roll,'r-')
hold on
plot(dt,pitch,'g-')
plot(dt,yaw,'b-')
plot([0 60],[pi/2 pi/2],'k-.')
plot([0 60],-[pi/2 pi/2],'k-.')
%{
plot(dt,roll-roll(end)/dt(end)*dt,'r-')
hold on
plot(dt,pitch-pitch(end)/dt(end)*dt,'g-')
plot(dt,yaw-yaw(end)/dt(end)*dt,'b-')
plot([0 60],[pi/4 pi/4],'k-.')
plot([0 60],-[pi/4 pi/4],'k-.')
%}
grid on
title('IMU Angle plot')
legend('roll','pitch','yaw')
ylabel('angle [rad]')
xlabel('time [s]')

params.ax_bias = ax_bias;
params.ay_bias = ay_bias;
params.az_bias = az_bias;
params.roll_bias = roll_bias;
params.pitch_bias = pitch_bias;
params.yaw_bias = yaw_bias;
params.sf_a = sf_a;
params.sf_w = sf_w;

save('imu_params.mat','params')