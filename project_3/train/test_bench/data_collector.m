clear all
close all
clc

%% orientation plotting initialization
plotting = false;
if plotting
    figure(1)
    clf
    x_plot = plot3(0,0,0,'r-');
    hold on
    y_plot = plot3(0,0,0,'g-');
    z_plot = plot3(0,0,0,'b-');
    axis equal
    grid on
    lims = [-1.1 1.1];
    xlim(lims)
    ylim(lims)
    zlim(lims)
end

%% load data
pattern = {'beat3','beat4','circle','eight','inf','wave'};
for pattern_num = 1:length(pattern)
    addpath(pattern{pattern_num})

    contents = dir(['../' pattern{pattern_num}]);
    file_names = cell(length(contents),1);
    for file_idx = 1:length(contents)
        file_names{file_idx} = contents(file_idx).name;
    end
    valid_files = find(~cellfun(@isempty,regexp(file_names,'.+\.txt')));
    for file_idx = 1:length(valid_files)
        file_name = file_names{valid_files(file_idx)};
        fid = fopen(['../' pattern{pattern_num} '/' file_name]);
        data = textscan(fid,'%d %f %f %f %f %f %f',...
            'TreatAsEmpty',{'NA','na'},'CommentStyle','#');
        fclose(fid);

        time = double(cell2mat(data(:,1)))/1000;
        acc = cell2mat(data(:,2:4));
        gyro = cell2mat(data(:,5:7));

        %% initialize
        % initial state mean and covariance
        x_k = [[1 0 0 0]'; zeros(3,1)];
        P_k = eye(6);

        % process and sensor noise
        Q = diag([78.98;78.98;78.98;18.9;18.9;18.9]);
        R = diag([0.5;0.5;0.5;0.001;0.001;0.001]);
        
        tic
        start = toc;
        O = [acc'; gyro'];
        
        %{
        O = zeros(4,length(time));
        O(:,1) = [1 0 0 0]';
        for i = 2:length(time)
            % get dt
            dt = time(i)-time(i-1);
            v_k = [(acc(i,:)/9.81)'; gyro(i,:)'];
            [x_k, P_k] = ukf(x_k,P_k,Q,R,v_k,dt);
            
            % get orientation
            O(:,i) = x_k(1:4);

            if plotting
                %% plot axes
                rot = quat_to_rot(x_k(1:4));
                coords = rot*eye(3);
                set(x_plot,'xdata',[coords(1,1) 0],'ydata',[coords(2,1) 0],'zdata',[coords(3,1) 0])
                set(y_plot,'xdata',[coords(1,2) 0],'ydata',[coords(2,2) 0],'zdata',[coords(3,2) 0])
                set(z_plot,'xdata',[coords(1,3) 0],'ydata',[coords(2,3) 0],'zdata',[coords(3,3) 0])
                stop = toc;
                pause(dt-(stop-start))
                start = stop;
            end
        end
        %}
        if ~plotting
            full_data.orientation = O;
            full_data.acceleration = acc';
            full_data.time = time;
            save(['data/' file_name(1:end-4) '.mat'],'full_data');
        end
    end
end