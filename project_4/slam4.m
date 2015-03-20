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

%% occupancy grid initialization
x_len = 30;
y_len = 30;
z_len = 5;
resolution = 0.05;
offset_x = x_len/2;
offset_y = y_len/2;
offset_z = 0;
dec_place = floor(log10(resolution));
val = floor(abs(resolution) ./ 10.^dec_place);
map_2 = ones(floor(x_len/resolution),floor(y_len/resolution));
path = ones(length(lidar.lidar),3)*0.5;

%% plot
figure(4)
clf
lidar_2d = imshow(1./(1+exp(-map_2))); % log odds map plot
hold on
grid on
set(gca,'ydir','normal')
axis equal
robot_pos = plot(offset_x,offset_y,'r-'); % robot path plot
robot_y = [1 0 -1]*0.5;
robot_x = [-1 2 -1]*0.5;
robot_plot = fill(robot_x,robot_y,'g'); % robot pose plot

%% time stamps of lidar data
lidar_t = zeros(1,length(lidar.lidar));
lidar_idx = 1;
for i = 1:length(lidar.lidar)
    lidar_t(i) = lidar.lidar{i}.t-joints.t0;
end
lidar_angles = linspace(-135,135,1081)'*pi/180;

%% main loop
nbytes = 0;
num_angles = 8;
patch_width = 3;
num_particles = num_angles*patch_width^2;
log_odd_threshold = 7;
H_bod = repmat(eye(4),1,1,num_particles);
lidar_points = zeros(1081,num_particles*4);
weights = ones(1,num_particles)/num_particles;

best_particle = zeros(3,1);
[o_x,o_y] = meshgrid(0:resolution:(patch_width-1)*resolution);
particle_offset = [reshape(o_x,1,patch_width^2)-floor(patch_width/2)*resolution;...
                    reshape(o_y,1,patch_width^2)-floor(patch_width/2)*resolution];
particles = zeros(3,num_particles);
[~,i_0] = min(abs(lidar_t(1) - joints.ts));
for i = i_0:length(joints.ts)
    %% process new lidar hits
    [~,idx] = min(abs(lidar_t - joints.ts(i)));
    if idx ~= lidar_idx
        lidar_idx = idx;

        %% get odometry data in local frame
        body_roll = lidar.lidar{lidar_idx}.rpy(1);
        body_pitch = lidar.lidar{lidar_idx}.rpy(2);
        body_yaw = lidar.lidar{lidar_idx}.rpy(3) - lidar.lidar{1}.rpy(3);
        bod_pose = lidar.lidar{lidar_idx}.pose;
        bod_pose(3) = body_yaw;
                
        if lidar_idx == 1
            dX = [0 0]';
            dtheta = 0;
        else
            body_yaw_prev = lidar.lidar{lidar_idx-1}.rpy(3) - lidar.lidar{1}.rpy(3);
            bod_pose_prev = lidar.lidar{lidar_idx-1}.pose;
            bod_pose_prev(3) = body_yaw_prev;
            dtheta = body_yaw-body_yaw_prev;
            dX = [cos(body_yaw_prev) sin(body_yaw_prev); -sin(body_yaw_prev) cos(body_yaw_prev)]*(bod_pose(1:2)'-bod_pose_prev(1:2)');
        end
        
        %% motion model update on particles
        thetas = dtheta + normrnd(0,abs(dtheta),[1 num_angles]);
        for j = 1:num_angles
            particles(1:2,(j-1)*patch_width^2+1:j*patch_width^2) = bsxfun(@plus,best_particle(1:2),particle_offset);
            particles(3,(j-1)*patch_width^2+1:j*patch_width^2) = best_particle(3) + thetas(j);
        end
        
        for j = 1:num_particles
            particles(1:2,j) = particles(1:2,j) + [cos(particles(3,j)) -sin(particles(3,j)); sin(particles(3,j)) cos(particles(3,j))]*dX + [normrnd(0,abs(dX(1)));normrnd(0,abs(dX(2)))];
        end
        
        %% get homogeneous transforms of each particle
        for j = 1:num_particles
            H_bod(:,:,j) = get_hom_transform(euler_to_rot(body_roll,body_pitch,particles(3,j)),[particles(1:2,j); 0]);
        end
        
        %% get homogeneous transform of head
        H_cam = get_hom_transform(eye(3),[0 0 0.395])*get_hom_transform(euler_to_rot(0,joints.pos(i,iHead),joints.pos(i,iNeck)),[0 0 0.085]);
        
        %% get rotated lidar scan in 3d (no translation)
        lidar_scan = lidar.lidar{lidar_idx}.scan';
        lidar_cart = [lidar_scan.*cos(lidar_angles) lidar_scan.*sin(lidar_angles) zeros(length(lidar_scan),1) ones(length(lidar_scan),1)];
        
        for j = 1:num_particles
            lidar_points(:,4*j-3:4*j) = (H_bod(:,:,j)*H_cam*lidar_cart')';
        end
        
        % get rid of bad lidar hits
        lidar_mask = lidar_scan>0.025;
        lidar_window = 1:1081;
        lidar_points_masked = lidar_points(lidar_mask(lidar_window),:);
        
        %% get correlations
        lidar_scaled = bsxfun(@plus,round(lidar_points_masked/resolution), repmat(round([offset_x offset_y offset_z 0]/resolution),1,num_particles));
        for j = 1:num_particles
            lidar_test = lidar_scaled(:,4*j-3:4*j-1);
            mask = sum(bsxfun(@gt, lidar_test(:,1:2), [x_len y_len]/resolution),2) | sum(bsxfun(@lt,lidar_test(:,1:2),[1 1]),2) | lidar_test(:,3) < 0.1/resolution;
            lidar_test(mask > 0,:) = [];
            lidar_linidx = sub2ind(size(map_2),lidar_test(:,1),lidar_test(:,2));
            weights(j) = sum(1./(1+exp(-map_2(lidar_linidx))));
        end
        weights = weights/sum(weights);
        
        %% update map based on most likely particle
        % get lidar hits of best particle
        [~,best_particle_idx] = max(weights);
        best_particle = particles(:,best_particle_idx);
        lidar_hits = lidar_scaled(:,best_particle_idx*4-3:best_particle_idx*4-2);
        lidar_z = lidar_scaled(:,best_particle_idx*4-1);
        
        % update positive hits
        mask = sum(bsxfun(@gt, lidar_hits, [x_len y_len]/resolution),2) | sum(bsxfun(@lt,lidar_hits,[1 1]),2) | lidar_z < 0.1/resolution | lidar_z > 1.5/resolution;
        lidar_hits(mask>0,:) = [];
        lidar_linidx = sub2ind(size(map_2),lidar_hits(:,1),lidar_hits(:,2));
        map_2(uint64(lidar_linidx)) = map_2(uint64(lidar_linidx)) + 1;
        
        % update negative hits
        x_start = round((best_particle(1)+offset_x)/resolution);
        y_start = round((best_particle(2)+offset_y)/resolution);
        [x_b,y_b] = getMapCellsFromRay(repmat(x_start,1,size(lidar_hits,1))',...
                                       repmat(y_start,1,size(lidar_hits,1))',...
                                        double(lidar_hits(:,1))',...
                                        double(lidar_hits(:,2))');
                                    
        mask = x_b > 0 & y_b > 0 & x_b < size(map_2,1) & y_b < size(map_2,2);
        lidar_nlinidx = sub2ind(size(map_2),x_b(mask),y_b(mask));
        map_2(uint64(lidar_nlinidx)) = map_2(uint64(lidar_nlinidx)) - 0.1;
        map_2(map_2 > log_odd_threshold) = log_odd_threshold;
        map_2(map_2 < -log_odd_threshold) = -log_odd_threshold;
        
        %% store path
        pos = [best_particle(1:2)' 0]';
        if sum(pos' == path(end,:)) < 3
            path(lidar_idx,:) = pos';
        end
    end
    
    %% redraw plots
    n = 20;
    if ~mod(i,n)
        % print the time
        fprintf(repmat('\b',1,nbytes));
        nbytes = fprintf('time: %6.6f',joints.ts(i));
        
        % plot 2d
        set(lidar_2d,'cdata',imrotate(1./(1+exp(-fliplr(map_2))),90));
        set(robot_pos,'xdata',(path(1:lidar_idx,1)+offset_x)/resolution,'ydata',(path(1:lidar_idx,2)+offset_y)/resolution)
        new_robot_pos = [cos(lidar.lidar{lidar_idx}.rpy(3) - lidar.lidar{1}.rpy(3)) -sin(lidar.lidar{lidar_idx}.rpy(3) - lidar.lidar{1}.rpy(3)); sin(lidar.lidar{lidar_idx}.rpy(3) - lidar.lidar{1}.rpy(3)) cos(lidar.lidar{lidar_idx}.rpy(3) - lidar.lidar{1}.rpy(3))]*[robot_x;robot_y];
        set(robot_plot,'xdata',(new_robot_pos(1,:)+path(lidar_idx,1)+offset_x)/resolution,'ydata',(new_robot_pos(2,:)+path(lidar_idx,2)+offset_y)/resolution)
        drawnow
    end
end
fprintf('\n')