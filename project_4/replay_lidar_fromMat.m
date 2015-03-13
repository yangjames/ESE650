load lidar0.mat

theta = 0:0.25:270;
theta = theta*pi/180;

for i=10:10:numel(lidar)
    % remove noisy data out of valid range
    lidar{i}.scan(find(lidar{i}.scan > 30)) = 0;
    
    figure(1), hold off;
    polar(theta, lidar{i}.scan,'b.');
    
    pause(0.025);
    i
end


