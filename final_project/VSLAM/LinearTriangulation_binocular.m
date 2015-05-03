function X = LinearTriangulation_binocular(K1,K2,C1,R1,C2,R2,x1,x2)
P1 = K1*R1*[eye(3) -C1];
P2 = K2*R2*[eye(3) -C2];
N = size(x1,1);

X = zeros(N,3);
for i = 1:N
    A = [cross_mat([x1(i,:) 1])*P1;...
        cross_mat([x2(i,:) 1])*P2];
    [~,~,V] = svd(A);
    X_h = V(:,end)/V(end,end);
    X(i,:) = X_h(1:3)';
end


    function C = cross_mat(v)
        C = [0 -v(3) v(2);...
            v(3) 0 -v(1);...
            -v(2) v(1) 0];
    end
end