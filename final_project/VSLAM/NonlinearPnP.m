function [C, R] = NonlinearPnP(X,x,K,C0,R0)

opts = optimoptions(@lsqnonlin, 'Algorithm', 'levenberg-marquardt',...
                                'TolX', 1e-64,...
                                'TolFun', 1e-64,...
                                'MaxFunEvals', 1e64,...
                                'MaxIter', 1e64,...
                                'Display','none');

[params,~] = lsqnonlin(@error,[R2q(R0);C0],[],[],opts,X,x,K);


R = q2R(params(1:4));
C = params(5:end);

    function err = error(param0, X, x, K)
        X_aug = [X ones(size(X,1),1)];
        P_n = K*q2R(param0(1:4))*[eye(3) -param0(5:end)];
        x_p = bsxfun(@rdivide,P_n(1:2,:)*X_aug',P_n(3,:)*X_aug')';
        err = x-x_p(:,1:2);
    end
end