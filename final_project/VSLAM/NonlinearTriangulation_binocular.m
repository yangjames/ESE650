function X = NonlinearTriangulation_binocular(K1, K2, C1, R1, C2, R2, x1, x2, X0)

opts = optimoptions(@lsqnonlin, 'Algorithm', 'levenberg-marquardt',...
                                'TolX', 1e-64,...
                                'TolFun', 1e-64,...
                                'MaxFunEvals', 1e64,...
                                'MaxIter', 1e64,...
                                'Display','none');

P1 = K1*R1*[eye(3) -C1];
P2 = K2*R2*[eye(3) -C2];
[N,M] = size(X0);
X = zeros([N,M]);
nbytes = 0;
for i = 1:size(X0,1)
    fprintf(repmat('\b',1,nbytes))
    nbytes = fprintf('iteration: %d of %d',i, N);
    [X(i,:),~] = lsqnonlin(@error,X0(i,:),[],[],opts,x1(i,:),x2(i,:),P1,P2);
end
fprintf('\n')
    function err = error(X1, x1, x2, P1, P2)
        X_aug = [X1 ones(size(X1,1),1)];
        x_p1 = bsxfun(@rdivide,(P1(1:2,:)*X_aug'),(P1(3,:)*X_aug'))';
        x_p2 = bsxfun(@rdivide,(P2(1:2,:)*X_aug'),(P2(3,:)*X_aug'))';
        err = [x1-x_p1 x2-x_p2];
    end
end