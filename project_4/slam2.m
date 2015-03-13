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

% image plot
f2 = figure(2);
clf
im = djpeg(rgb.RGB{1}.jpeg);
im_plot = imagesc(flip(im,2));

% grid plot
figure(3)
clf
x_att_plot = plot3(0,0,0,'r-');
hold on
grid on
axis equal
y_att_plot = plot3(0,0,0,'g-');
z_att_plot = plot3(0,0,0,'b-');
axes_scale = 1;
path_plot = plot3(0,0,0,'k-');
lidar_plot = plot3(0,0,0,'b.');

% map plot
figure(4)
clf
robot_y = [1 0 -1]*0.5;
robot_x = [-1 2 -1]*0.5;
robot_plot = fill(robot_x,robot_y,'g');
hold on
grid on
axis equal
robot_pos = plot(0,0,'r-');
lidar_2d = plot(0,0,'b.','MarkerSize',1);

% grid initialization
x_len = 20;
y_len = 20;
z_len = 3;
resolution = 0.05;
offset_x = x_len/2;
offset_y = y_len/2;
offset_z = z_len/2;

dec_place = floor(log10(resolution));
val = floor(abs(resolution) ./ 10.^dec_place);

occupancy_grid = false(floor(x_len/resolution),floor(y_len/resolution),floor(z_len/resolution));

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

% initial body attitude
R_bod = eye(3);
pos = [0 0 0];
path = [0 0 0];

% points inside a 10m circle from 0 to 270 degrees
r = 0:0.01:5;
in_points = zeros(length(r)*length(lidar_angles),2);
idx = 1;
for i = 1:length(r)
    for j = 1:length(lidar_angles)
        in_points(idx,:) = [r(i) lidar_angles(j)];
        idx = idx+1;
    end
end

% main loop
for i = 2:length(joints.ts)
    % draw local pose of head
    yaw_h = joints.pos(i,iNeck);
    pitch_h = joints.pos(i,iHead);
    R_cam = R_bod*euler_to_rot(0,pitch_h,yaw_h);
    set(x_att_plot,'xdata',[0 R_cam(1,1)]*axes_scale + pos(1),...
                    'ydata',[0 R_cam(2,1)]*axes_scale + pos(2),...
                    'zdata',[0 R_cam(3,1)]*axes_scale)
    set(y_att_plot,'xdata',[0 R_cam(1,2)]*axes_scale + pos(1),...
                    'ydata',[0 R_cam(2,2)]*axes_scale + pos(2),...
                    'zdata',[0 R_cam(3,2)]*axes_scale)
    set(z_att_plot,'xdata',[0 R_cam(1,3)]*axes_scale + pos(1),...
                    'ydata',[0 R_cam(2,3)]*axes_scale + pos(2),...
                    'zdata',[0 R_cam(3,3)]*axes_scale)
    
    % draw the lidar hits that we see
    [~,idx] = min(abs(lidar_t - joints.ts(i)));
    if idx ~= lidar_idx
        % get closest matching lidar index to time stamp
        lidar_idx = idx;

        % get rotation matrix from world to body frame
        R_bod = euler_to_rot(lidar.lidar{lidar_idx}.rpy(1),lidar.lidar{lidar_idx}.rpy(2),lidar.lidar{lidar_idx}.pose(3));
        
        % get xy position in world
        pos = [lidar.lidar{lidar_idx}.pose(1:2) 0];
        if sum(pos == path(end,:)) < 3
            path = [path;pos(1:2) 0];
            set(path_plot,'xdata',path(:,1),'ydata',path(:,2),'zdata',path(:,3))
        end
        
        % get lidar values
        lidar_scan = lidar.lidar{lidar_idx}.scan';
        
        % get lidar hits in cartesian space
        lidar_points = [[lidar_scan.*cos(lidar_angles) lidar_scan.*sin(lidar_angles) zeros(length(lidar_scan),1)]; [0 0 0]];
        lidar_points = lidar_points(lidar_scan>0.025,:);
        
        %{
        %get cells inside lidar contour
        lidar_in = inpolygon(lidar_points(:,1),lidar_points(:,2),in_points(:,1),in_points(:,2));
        %}
        
        % get the lidar hits in cartesian space wrt the world
        lidar_pose = lidar.lidar{lidar_idx}.pose;
        lidar_rot = bsxfun(@plus,round((R_cam*lidar_points')/val,int8(-dec_place))*val,[lidar_pose(1:2) 0]');
        %lidar_in_rot = bsxfun(@plus,round((R_cam*lidar_in')/val,int8(-dec_place))*val,[lidar_pose(1:2) 0]');
        
        % scale lidar hits for the map
        lidar_rot_scaled = bsxfun(@plus,round(lidar_rot/resolution),round([offset_x offset_y offset_z]'/resolution));
        %lidar_in_rot_scaled = bsxfun(@plus,round(lidar_rot/resolution),round([offset_x offset_y offset_z]'/resolution));
        
        % prune the lidar hits that exceed map dimensions
        mask = sum(bsxfun(@gt, lidar_rot_scaled, [x_len y_len z_len]'/resolution));
        lidar_rot_scaled(:,mask > 0) = [];
        mask = sum(bsxfun(@lt, lidar_rot_scaled, [1 1 1]'));
        lidar_rot_scaled(:,mask > 0) = [];
        
        %{
        mask = sum(bsxfun(@gt, lidar_in_rot_scaled, [x_len y_len z_len]'/resolution));
        lidar_in_rot_scaled(:,mask > 0) = [];
        mask = sum(bsxfun(@lt, lidar_in_rot_scaled, [1 1 1]'));
        lidar_in_rot_scaled(:,mask > 0) = [];
        %}
        
        % fill grid with new hits
        lidar_linidx = sub2ind(size(occupancy_grid),lidar_rot_scaled(1,:),lidar_rot_scaled(2,:),lidar_rot_scaled(3,:));
        occupancy_grid(uint64(lidar_linidx)) = true;
        
        %{
        % empty grid cells that clearly cannot be filled (within line of
        % sight of the lidar
        not_lidar_linidx = sub2ind(size(occupancy_grid),lidar_in_rot_scaled(1,:),lidar_rot_scaled(2,:),lidar_rot_scaled(3,:));
        occupancy_grid(uint64(lidar_linidx)) = false;
        %}
        
        % plot updated grid
        populated = find(occupancy_grid);
        [x,y,z] = ind2sub(size(occupancy_grid),populated);
        set(lidar_plot,'xdata',x*resolution-offset_x,'ydata',y*resolution-offset_y,'zdata',z*resolution-offset_z)
        
        % plot updated 2d map
        set(lidar_2d,'xdata',x*resolution-offset_x,'ydata',y*resolution-offset_y)
        set(robot_pos,'xdata',path(:,1),'ydata',path(:,2))
        new_robot_pos = [cos(lidar.lidar{lidar_idx}.pose(3)) -sin(lidar.lidar{lidar_idx}.pose(3)); sin(lidar.lidar{lidar_idx}.pose(3)) cos(lidar.lidar{lidar_idx}.pose(3))]*[robot_x;robot_y];
        set(robot_plot,'xdata',new_robot_pos(1,:)+pos(1),'ydata',new_robot_pos(2,:)+pos(2))
    end
    [~,idx] = min(abs(rgb_t-joints.ts(i)));
    disp([rgb_t(idx) joints.ts(i)])
    
    % draw every nth frame
    n = 10;
    if ~mod(i,n)
        drawnow
    end
end