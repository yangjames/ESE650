function [x_k,P_k] = ukf(x_k,P_k,Q,R,v_k,dt)
    X = zeros(7,13);
    Y = X;

    n = length(P_k);
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
    v_k = v_k-z_k_mean;
    v_kc = K_k*v_k;

    alpha_v = norm(v_kc(1:3));
    e_v = v_kc(1:3)/norm(v_kc(1:3));
    e_v(:,alpha_v == 0) = repmat([0 0 0]',1,sum(alpha_v==0));
    q_v = [cos(alpha_v/2);...
        e_v*sin(alpha_v/2)];
    q_combined = bsxfun(@quat_mult, x_k_mean(1:4),q_v);
    x_k = [bsxfun(@rdivide,q_combined,sqrt(sum(q_combined.^2))); x_k_mean(5:7)+v_kc(4:6)];
    P_k = P_k_mean-K_k*P_vv*K_k';

    function prod = quat_mult(q,p)
        prod = [q(1)*p(1)-q(2)*p(2)-q(3)*p(3)-q(4)*p(4);...
        q(2)*p(1)+q(1)*p(2)-q(4)*p(3)+q(3)*p(4);...
        q(3)*p(1)+q(4)*p(2)+q(1)*p(3)-q(2)*p(4);...
        q(4)*p(1)-q(3)*p(2)+q(2)*p(3)+q(1)*p(4)];
    end
    function q_bar = quat_conj(q)
        q_bar = bsxfun(@times,q,[1;-1;-1;-1]);
    end
end