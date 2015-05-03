function [y1, y2, idx, H_best] = GetPlaneRANSAC(X1_g,x2_g,threshold, max_iter)

N = size(X1_g,1);
y1 = [];
y2 = [];
idx = [];
H_best = [];
if N < 4
    return
end
n_best = 0;
for i = 1:max_iter
    r_idx = ceil(rand(4,1)*N);
    H=EstimateHomography(x2_g(r_idx,1),x2_g(r_idx,2), X1_g(r_idx,1),X1_g(r_idx,2));
    trans = H*X1_g';
    mask = sqrt(sum((bsxfun(@rdivide,trans(1:2,:),trans(3,:))-x2_g').^2)) < threshold;
    num_in = sum(mask);
    if n_best < num_in
        n_best = num_in;
        H_best = H;
        y1 = X1_g(mask,:);
        y2 = x2_g(mask,:);
        idx = mask;%find(mask);
    end
end