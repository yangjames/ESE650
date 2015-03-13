clear all
close all

% add important directories
addpath('mex')
addpath('cameraParam')

% load data
joints = load('joints0.mat');
lidar = load('lidar0.mat');
rgb = load('rgb0.mat');
depth = load('depth0.mat');
idx = 1;
iNeck = get_joint_index('Neck'); % head yaw
iHead = get_joint_index('Head'); % head pitch

% initialize plotting
figure(2)
plot(joints.ts, joints.pos(:, iNeck));

figure(3)
clf
x_att_plot = plot3(0,0,0,'r-');
hold on
grid on
axis equal
y_att_plot = plot3(0,0,0,'g-');
z_att_plot = plot3(0,0,0,'b-');
path_plot = plot3(0,0,0,'k-');
lidar_plot = plot3(0,0,0,'b.');

% get first image
f2 = figure(2);
clf
im = djpeg(rgb.RGB{1}.jpeg);
im_plot = imagesc(flip(im,2));

% get camera timestamps
rgb_t = [];
im_idx = 1;
for i = 1:length(rgb.RGB)
    rgb_t = [rgb_t rgb.RGB{i}.t-joints.t0];
end

% get lidar timestamps
lidar_t = [];
lidar_idx = 1;
for i = 1:length(lidar.lidar)
    lidar_t = [lidar_t lidar.lidar{i}.t-joints.t0];
end
lidar_angles = linspace(-135,135,1081)'*pi/180;

% initialize storage
pos = [0 0 0];
path = [0 0 0];
R_bod = eye(3);
resolution = 0.05;
dec_place = floor(log10(resolution));
val = floor(abs(resolution) ./ 10.^dec_place);
x_len = 20;
y_len = 20;
z_len = 3;
offset_x = x_len/2;
offset_y = y_len/2;
offset_z = z_len/2;

occupancy_grid = false(floor(x_len/resolution),floor(y_len/resolution),floor(z_len/resolution));
scale = 1;
max_x = 1;
max_y = 1;
max_z = 1;
min_x = 0;
min_y = 0;
min_z = 0;
% main loop
for i = 2:length(joints.ts)
    % draw local pose of head
    yaw_h = joints.pos(i,iNeck);
    pitch_h = joints.pos(i,iHead);
    R_cam = R_bod*euler_to_rot(0,pitch_h,yaw_h);
    set(x_att_plot,'xdata',[0 R_cam(1,1)]*scale + pos(1),'ydata',[0 R_cam(2,1)]*scale + pos(2),'zdata',[0 R_cam(3,1)]*scale)
    set(y_att_plot,'xdata',[0 R_cam(1,2)]*scale + pos(1),'ydata',[0 R_cam(2,2)]*scale + pos(2),'zdata',[0 R_cam(3,2)]*scale)
    set(z_att_plot,'xdata',[0 R_cam(1,3)]*scale + pos(1),'ydata',[0 R_cam(2,3)]*scale + pos(2),'zdata',[0 R_cam(3,3)]*scale)

    % draw the lidar hits that we see
    [~,idx] = min(abs(lidar_t - joints.ts(i)));
    if idx ~= lidar_idx
        lidar_idx = idx;
        lidar_scan = lidar.lidar{lidar_idx}.scan';
        lidar_points = [lidar_scan.*cos(lidar_angles) lidar_scan.*sin(lidar_angles) zeros(length(lidar_scan),1)];
        lidar_pose = lidar.lidar{lidar_idx}.pose;
        lidar_rot = bsxfun(@plus,round((R_cam*lidar_points')/val,int8(-dec_place))*val,[lidar_pose(1:2) 0]');
        
        lidar_rot_scaled = bsxfun(@plus,round(lidar_rot/resolution),round([offset_x offset_y offset_z]'/resolution));
        
        mask = sum(bsxfun(@gt, lidar_rot_scaled, [x_len y_len z_len]'/resolution));
        lidar_rot_scaled(:,mask > 0) = [];
        mask = sum(bsxfun(@lt, lidar_rot_scaled, [1 1 1]'));
        lidar_rot_scaled(:,mask > 0) = [];
        
        lidar_linidx = sub2ind(size(occupancy_grid),lidar_rot_scaled(1,:),lidar_rot_scaled(2,:),lidar_rot_scaled(3,:));
        occupancy_grid(uint64(lidar_linidx)) = true;
        
        populated = find(occupancy_grid);
        [x,y,z] = ind2sub(size(occupancy_grid),populated);
        
        %{
        if max_x < max(x*resolution-offset_x)
            max_x = max(x*resolution-offset_x);
        end
        if max_y < max(y*resolution-offset_y)
            max_y = max(y*resolution-offset_y);
        end
        if max_z < max(z*resolution-offset_z)
            max_z = max(z*resolution-offset_z);
        end
        
        if min_x > min(x*resolution-offset_x)
            min_x = min(x*resolution-offset_x);
        end
        if min_y > min(y*resolution-offset_y)
            min_y = min(y*resolution-offset_y);
        end
        if max_z > min(z*resolution-offset_z)
            min_z = min(z*resolution-offset_z);
        end
        %set(lidar_plot,'xlim',[min_x max_x])
        xlim([min_x max_x])
        ylim([min_y max_y])
        zlim([min_z max_z])
        %}
        set(lidar_plot,'xdata',x*resolution-offset_x,'ydata',y*resolution-offset_y,'zdata',z*resolution-offset_z)
        %set(lidar_plot,'xdata',lidar_rot_scaled(1,:)*resolution-offset_x,'ydata',lidar_rot_scaled(2,:)*resolution-offset_y,'zdata',lidar_rot_scaled(3,:)*resolution-offset_z)
    end
    
    % draw the image that we see
    [~,idx] = min(abs(rgb_t-joints.ts(i)));
    %disp([rgb_t(idx) joints.ts(i)])
    if idx ~= im_idx
        im_idx = idx;
        im = djpeg(rgb.RGB{im_idx}.jpeg);
        pos = rgb.RGB{im_idx}.odom;
        if sum(pos == path(end,:)) < 3
            path = [path;pos(1:2) 0];
            set(path_plot,'xdata',path(:,1),'ydata',path(:,2),'zdata',path(:,3))
        end
        
        imu_rpy = rgb.RGB{im_idx}.imu_rpy;
        R_bod = euler_to_rot(imu_rpy(1),imu_rpy(2),lidar_pose(3));
        set(im_plot,'cdata',flip(im,2))
    end
    
    % draw every fifth frame
    if mod(i,10) == 0
        drawnow
    end
end