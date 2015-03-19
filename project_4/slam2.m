clear all
close all

%% add important directories
addpath('mex')
addpath('cameraParam')

%% load data
joints = load('joints0.mat');
lidar = load('lidar0.mat');
rgb = load('rgb0.mat');
depth = load('depth0.mat');
iNeck = get_joint_index('Neck'); % head yaw
iHead = get_joint_index('Head'); % head pitch

%% plot flags
plot3d = false;
show_images = false;

%% grid initialization
x_len = 40;
y_len = 40;
z_len = 5;
resolution = 0.1;
offset_x = x_len/2;
offset_y = y_len/2;
offset_z = z_len/2;

dec_place = floor(log10(resolution));
val = floor(abs(resolution) ./ 10.^dec_place);

if plot3d
    occupancy_grid_3 = false(floor(x_len/resolution),floor(y_len/resolution),floor(z_len/resolution));
end

map_2 = ones(floor(x_len/resolution),floor(y_len/resolution));

%% plots
% image plot
if show_images
    f2 = figure(2);
    clf
    im = djpeg(rgb.RGB{1}.jpeg);
    im_plot = imagesc(flip(im,2));
    
    % get camera timestamps
    rgb_t = zeros(1,length(rgb.RGB));
    im_idx = 1;
    for i = 1:length(rgb.RGB)
        rgb_t(i) = rgb.RGB{i}.t-joints.t0;
    end
end

% grid plot
if plot3d
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
    lidar_plot = plot3(0,0,0,'b.','MarkerSize',0.5);
end

% map plot
figure(4)
clf
lidar_2d = imshow(1./(1+exp(-map_2)));
hold on
grid on
set(gca,'ydir','normal')
axis equal
robot_y = [1 0 -1]*0.5;
robot_x = [-1 2 -1]*0.5;
robot_plot = fill(robot_x,robot_y,'g');
robot_pos = plot(offset_x,offset_y,'r-');
set(gca,'XTickLabel',get(gca,'XTick')*resolution)

%% time stamps of pictures and lidar data
% get lidar timestamps
lidar_t = zeros(1,length(lidar.lidar));
lidar_idx = 1;
for i = 1:length(lidar.lidar)
    lidar_t(i) = lidar.lidar{i}.t-joints.t0;
end
lidar_angles = linspace(-135,135,1081)'*pi/180;

% initial body attitude
R_bod = eye(3);
pos = [0 0 0];
path = ones(length(lidar.lidar),3)*0.5;

% main loop
nbytes = 0;
p_lidar = 0.55;
num_particles = 50;

for i = 2:length(joints.ts)

    % draw the lidar hits that we see
    [~,idx] = min(abs(lidar_t - joints.ts(i)));
    if idx ~= lidar_idx
        % get closest matching lidar index to time stamp
        lidar_idx = idx;

        % get previous motion
        bod_rpy = lidar.lidar{lidar_idx}.rpy - [0 0 lidar.lidar{1}.rpy(3)];
        bod_pose = lidar.lidar{lidar_idx}.pose;
        bod_pose(3) = bod_rpy(3);
        
        if lidar_idx - 1 > 0
            bod_rpy_prev = lidar.lidar{lidar_idx-1}.rpy - [0 0 lidar.lidar{1}.rpy(3)];
            bod_pose_prev = lidar.lidar{lidar_idx-1}.pose;
            bod_pose_prev(3) = bod_rpy_prev(3);
        else
            bod_rpy_prev = zeros(1,3);
            bod_pose_prev = zeros(1,3);
        end
        R_bod_prev = euler_to_rot(bod_rpy_prev(1), bod_rpy_prev(2),bod_rpy_prev(3));
        dx = bod_pose(1)-bod_pose_prev(1);
        dy = bod_pose(2)-bod_pose_prev(2);
        dtheta = bod_pose(3)-bod_pose_prev(3);
        
        dpos = euler_to_rot(0,0,bod_rpy_prev(3))'*[dx; dy; 0];
        
        % get rotation matrix from world to body frame
        R_bod = euler_to_rot(bod_rpy(1),bod_rpy(2),bod_rpy(3));
        
        % get local pose of head
        yaw_h = joints.pos(i,iNeck);
        pitch_h = joints.pos(i,iHead);
        R_cam = R_bod*euler_to_rot(0,pitch_h,yaw_h);
        
        % get xy position in world
        pos = [lidar.lidar{lidar_idx}.pose(1:2) 0]';
        if sum(pos' == path(end,:)) < 3
            path(lidar_idx,:) = pos';
        end
        
        % get lidar values
        lidar_scan = lidar.lidar{lidar_idx}.scan';
        
        % get lidar hits in cartesian space
        lidar_points = [[lidar_scan.*cos(lidar_angles) lidar_scan.*sin(lidar_angles) zeros(length(lidar_scan),1)]; [0 0 0]];
        lidar_points = lidar_points(lidar_scan>0.025,:);
                        
        % get the lidar hits in cartesian space wrt the world
        lidar_pose = lidar.lidar{lidar_idx}.pose;
        lidar_rot = bsxfun(@plus,round((R_cam*lidar_points')/val,int8(-dec_place))*val,[lidar_pose(1:2) 0]');
                
        % scale lidar hits for the map
        lidar_rot_scaled = bsxfun(@plus,round((lidar_rot)/resolution),round([offset_x offset_y offset_z]'/resolution));
                
        % prune the lidar hits that exceed map dimensions
        mask = sum(bsxfun(@gt, lidar_rot_scaled, [x_len y_len z_len]'/resolution));
        lidar_rot_scaled(:,mask > 0) = [];
        mask = sum(bsxfun(@lt, lidar_rot_scaled, [1 1 1]'));
        lidar_rot_scaled(:,mask > 0) = [];
        
        if plot3d
            % fill grid with new hits
            lidar_linidx_3 = sub2ind(size(occupancy_grid_3),lidar_rot_scaled(1,:),lidar_rot_scaled(2,:),lidar_rot_scaled(3,:));
            occupancy_grid_3(uint64(lidar_linidx_3)) = true;
        end
        
        lidar_linidx_2 = sub2ind(size(map_2),lidar_rot_scaled(1,:),lidar_rot_scaled(2,:));
        map_2(uint64(lidar_linidx_2)) = map_2(uint64(lidar_linidx_2)) + log(p_lidar/(1-p_lidar));
        
        % empty grid cells that go past the lidar hits
        x_start = round((lidar_pose(1)+offset_x)/resolution);
        y_start = round((lidar_pose(2)+offset_y)/resolution);
        [x_b,y_b] = getMapCellsFromRay(repmat(x_start,1,size(lidar_rot_scaled,2)),...
                                       repmat(y_start,1,size(lidar_rot_scaled,2)),...
                                        double(lidar_rot_scaled(1,:)),...
                                        double(lidar_rot_scaled(2,:)));
        mask = x_b > 0 & y_b > 0 & x_b < size(map_2,1) & y_b < size(map_2,2);
        lidar_nlinidx = sub2ind(size(map_2),x_b(mask),y_b(mask));
        map_2(uint64(lidar_nlinidx)) = map_2(uint64(lidar_nlinidx)) - log(p_lidar/(1-p_lidar));
        map_2(map_2 > 1000) = 1000;
        map_2(map_2 < -1000) = -1000;
    end
    
    % draw every nth frame
    n = 20;
    if ~mod(i,n)
        % print the time
        fprintf(repmat('\b',1,nbytes));
        nbytes = fprintf('time: %6.6f',joints.ts(i));
        
        % display camera image
        if show_images
            % find the time stamp for the camera
            [~,im_idx] = min(abs(rgb_t-joints.ts(i)));
            im = djpeg(rgb.RGB{im_idx}.jpeg);
            set(im_plot,'cdata',flip(im,2))
        end
        
        % plot 3d
        if plot3d
            populated = find(occupancy_grid_3);
            [x,y,z] = ind2sub(size(occupancy_grid_3),populated);
            set(lidar_plot,'xdata',x*resolution-offset_x,'ydata',y*resolution-offset_y,'zdata',z*resolution-offset_z)

            set(x_att_plot,'xdata',[0 R_cam(1,1)]*axes_scale + pos(1),...
                            'ydata',[0 R_cam(2,1)]*axes_scale + pos(2),...
                            'zdata',[0 R_cam(3,1)]*axes_scale)
            set(y_att_plot,'xdata',[0 R_cam(1,2)]*axes_scale + pos(1),...
                            'ydata',[0 R_cam(2,2)]*axes_scale + pos(2),...
                            'zdata',[0 R_cam(3,2)]*axes_scale)
            set(z_att_plot,'xdata',[0 R_cam(1,3)]*axes_scale + pos(1),...
                            'ydata',[0 R_cam(2,3)]*axes_scale + pos(2),...
                            'zdata',[0 R_cam(3,3)]*axes_scale)
            set(path_plot,'xdata',path(:,1),'ydata',path(:,2),'zdata',path(:,3))
        end
        
        % plot 2d
        set(lidar_2d,'cdata',imrotate(1./(1+exp(-fliplr(map_2))),90));
        set(robot_pos,'xdata',(path(1:lidar_idx,1)+offset_x)/resolution,'ydata',(path(1:lidar_idx,2)+offset_y)/resolution)
        new_robot_pos = [cos(lidar.lidar{lidar_idx}.rpy(3) - lidar.lidar{1}.rpy(3)) -sin(lidar.lidar{lidar_idx}.rpy(3) - lidar.lidar{1}.rpy(3)); sin(lidar.lidar{lidar_idx}.rpy(3) - lidar.lidar{1}.rpy(3)) cos(lidar.lidar{lidar_idx}.rpy(3) - lidar.lidar{1}.rpy(3))]*[robot_x;robot_y];
        set(robot_plot,'xdata',(new_robot_pos(1,:)+path(lidar_idx,1)+offset_x)/resolution,'ydata',(new_robot_pos(2,:)+path(lidar_idx,2)+offset_y)/resolution)
        drawnow
    end
end
fprintf('\n')