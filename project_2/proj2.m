clear all
addpath('imu')
addpath('cam')
addpath('ref')
addpath('vicon')

dataset = 1;
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
az_bias = 500;
sf_a = 9.81/106;
figure(2)
clf
plot(imu.ts-imu.ts(1),(imu.vals(1,:)-ax_bias)*sf_a,'r-')
hold on
plot(imu.ts-imu.ts(1),(imu.vals(2,:)-ay_bias)*sf_a,'g-')
plot(imu.ts-imu.ts(1),-(imu.vals(3,:)-az_bias)*sf_a,'b-')
grid on
title('IMU acceleration plot')
legend('x','y','z')
ylabel('acceleration [m/s^2]')
xlabel('time [s]')

% gyro data
roll_bias = 373;
pitch_bias = 375;
yaw_bias = 369;
sf_w = 3300/1023*pi/180/3.7;
figure(3)
clf
plot(imu.ts-imu.ts(1),(imu.vals(5,:)-roll_bias)*sf_w,'r-')
hold on
plot(imu.ts-imu.ts(1),(imu.vals(6,:)-pitch_bias)*sf_w,'g-')
plot(imu.ts-imu.ts(1),(imu.vals(4,:)-yaw_bias)*sf_w,'b-')
grid on
title('IMU gyro plot')
legend('roll','pitch','yaw')
ylabel('angular rate [rad/s]')
xlabel('time [s]')

% angular rate integration plot
dt = imu.ts(2:end)-imu.ts(1:end-1);
roll = cumsum((imu.vals(5,2:end)-roll_bias)*sf_w.*dt);
pitch = cumsum((imu.vals(6,2:end)-pitch_bias)*sf_w.*dt);
yaw = cumsum((imu.vals(4,2:end)-yaw_bias)*sf_w.*dt);
figure(4)
clf
plot(imu.ts(2:end)-imu.ts(1),roll,'r-')
hold on
plot(imu.ts(2:end)-imu.ts(1),pitch,'g-')
plot(imu.ts(2:end)-imu.ts(1),yaw,'b-')
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