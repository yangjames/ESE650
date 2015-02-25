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

    contents = dir(pattern{pattern_num});
    file_names = cell(length(contents),1);
    for file_idx = 1:length(contents)
        file_names{file_idx} = contents(file_idx).name;
    end
    valid_files = find(~cellfun(@isempty,regexp(file_names,'.+\.txt')));
    for file_idx = 1:length(valid_files)
        file_name = file_names{valid_files(file_idx)};
        fid = fopen(file_name);
        data = textscan(fid,'%d %f %f %f %f %f %f',...
            'TreatAsEmpty',{'NA','na'},'CommentStyle','#');
        fclose(fid);

        time = double(cell2mat(data(:,1)))/1000;
        acc = cell2mat(data(:,2:4))/9.81;
        gyro = cell2mat(data(:,5:7));

        %% initialize
        % initial state mean and covariance
        x_k = [[1 0 0 0]'; zeros(3,1)];
        P_k = eye(6);
        n = length(P_k);

        % process and sensor noise
        Q = diag([78.98;78.98;78.98;18.9;18.9;18.9]);
        R = diag([0.5;0.5;0.5;0.001;0.001;0.001]);

        % intialize storage vectors
        Y = zeros(7,13);
        X = Y;
        Z = zeros(6,13);

        tic
        start = toc;
        orientation = zeros(4,length(time));
        orientation(:,1) = [1 0 0 0]';
        
        for i = 2:length(time)
            % get dt
            dt = time(i)-time(i-1);
            q_k = x_k(1:4);
            omega_k = x_k(end-2:end);

            % obtain sigma points using Cholesky decomposition
            S = chol(P_k+Q);
            W = [sqrt(n)*S zeros(6,1) -sqrt(n)*S];

            % propagate sigma points
            alpha_W = sqrt(sum(W(1:3,:).^2));
            e_W = bsxfun(@rdivide,W(1:3,:),alpha_W);
            e_W(:,alpha_W == 0) = repmat([0 0 0]',1,sum(alpha_W==0));
            q_W = [cos(alpha_W/2);...
                bsxfun(@times,e_W,sin(alpha_W/2))];

            for j = 1:size(q_W,2)
                X(:,j) = [quat_mult(q_k,q_W(:,j)); omega_k+W(4:6,j)];
            end
            X(1:4,:) = bsxfun(@rdivide,X(1:4,:),sqrt(sum(X(1:4,:).^2)));

            % transform sigma points
            alpha_delta = sqrt(sum(X(5:7,:).^2))*dt;
            e_delta = bsxfun(@rdivide,X(5:7,:)*dt,alpha_delta);
            e_delta(:,alpha_delta == 0) = repmat([0 0 0]',1,sum(alpha_delta==0));
            q_delta = [cos(alpha_delta/2);...
                bsxfun(@times,e_delta,sin(alpha_delta/2))];
            q_k1 = bsxfun(@quat_mult,q_k,q_delta);

            for j = 1:size(q_delta,2)
                Y(:,j) = [quat_mult(X(1:4,j),q_delta(:,j)); X(5:7,j)];
            end

            Y(1:4,:) = bsxfun(@rdivide,Y(1:4,:),sqrt(sum(Y(1:4,:).^2)));

            % calculate sigma point mean and covariance
            [~,~,V] = svd((Y(1:4,:)*Y(1:4,:)')/(2*n));
            x_k_mean = [V(:,1)/norm(V(:,1));mean(Y(end-2:end,:),2)];

            r_prime = bsxfun(@quat_mult,quat_conj(x_k_mean(1:4)),Y(1:4,:));
            omega_prime = bsxfun(@minus,Y(5:7,:),x_k_mean(5:7));
            Y_mean_centered = [r_prime; omega_prime];
            W_prime = [bsxfun(@times,2*acos(Y_mean_centered(1,:)) ...
                                    ./sin(acos(Y_mean_centered(1,:))),Y_mean_centered(2:4,:));...
                    Y_mean_centered(end-2:end,:)];
            mask = sin(acos(Y_mean_centered(1,:))) == 0;
            W_prime(1:4,mask) = repmat([1 0 0 0]',1,sum(mask));
            P_k_mean = W_prime*W_prime'/(2*n);

            % compute a-priori estimate
            g_rot = bsxfun(@quat_mult,quat_conj(Y(1:4,:)),[0 0 0 1]');
            for j = 1:size(Y,2)
                g_rot(:,j) = quat_mult(g_rot(:,j),Y(1:4,j));
            end

            Z = [g_rot(2:4,:);Y(end-2:end,:)];
            z_k_mean = mean(Z,2);

            Z_mean_centered = bsxfun(@minus,Z,z_k_mean);
            P_zz = (Z_mean_centered*Z_mean_centered')/(2*n);
            P_vv = P_zz + R;
            P_xz = (W_prime*Z_mean_centered')/(2*n);

            % calculate Kalman gain, residual, and state update
            K_k = P_xz/P_vv;
            v_k = [acc(i,:)';...
                gyro(i,:)']-z_k_mean;
            v_kc = K_k*v_k;

            x_q = quat_to_vec(x_k_mean(1:4)');

            alpha_v = norm(v_kc(1:3));
            e_v = v_kc(1:3)/norm(v_kc(1:3));
            e_v(:,alpha_v == 0) = repmat([0 0 0]',1,sum(alpha_v==0));
            q_v = [cos(alpha_v/2);...
                e_v*sin(alpha_v/2)];
            q_combined = bsxfun(@quat_mult, x_k_mean(1:4),q_v);
            x_k = [bsxfun(@rdivide,q_combined,sqrt(sum(q_combined.^2))); x_k_mean(5:7)+v_kc(4:6)];
            P_k = P_k_mean-K_k*P_vv*K_k';

            %% store orientation
            orientation(:,i) = x_k(1:4);

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
        if ~plotting
            full_data.orientation = orientation;
            full_data.acceleration = acc';
            full_data.time = time;
            save([file_name(1:end-4) '.mat'],'full_data');
        end
    end
end