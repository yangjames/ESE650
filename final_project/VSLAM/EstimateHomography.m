function H = EstimateHomography(x_to,y_to,X_from,Y_from)
A = zeros(length(X_from(:))*2,9);

for i = 1:length(X_from(:)),
 a = [X_from(i),Y_from(i),1];
 b = [0 0 0];
 c = [x_to(i);y_to(i)];
 d = -c*a;
 A((i-1)*2+1:(i-1)*2+2,1:9) = [[a b;b a] d];
end

[~, ~, V] = svd(A);
h = V(:,9);
H = reshape(h,3,3)';

H = H./H(3, 3);