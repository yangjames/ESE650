function [C, R] = LinearPnP(X, x, K)
N = length(X);

M = zeros(3*N,12);
for i=1:N
    M(3*(i-1)+1:3*i,:) = cross_mat(X(i,:),x(i,:));

end

[~,~,V] = svd(M);
P = reshape(V(:,end)/V(end,end),4,3)';
H = K\P;

[UR,~,UV] = svd(H(:,1:3));

R = sign(det(UR*UV'))*UR*UV';
C = -R'*H(:,4);

    function A = cross_mat(X_i,x_i)
        x_p = [x_i(:); 1];
        X_p = [X_i(:); 1];
        A = [0 -x_p(3) x_p(2);...
            x_p(3) 0 -x_p(1);...
            -x_p(2) x_p(1) 0]*[X_p' zeros(1,8);...
                            zeros(1,4) X_p' zeros(1,4);...
                            zeros(1,8) X_p'];
    end
end