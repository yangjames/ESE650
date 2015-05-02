function [y1,y2,idx] = GetInliersRANSAC(x1,x2,threshold, max_iter)
N = size(x1,1);
y1 = [];
y2 = [];
idx = [];
if N < 8
    return
end
x1_aug = [x1 ones(N,1)];
x2_aug = [x2 ones(N,1)];

n_best = 0;
for i = 1:max_iter
    r_idx = ceil(rand(8,1)*N);
    F=EstimateFundamentalMatrix(x1(r_idx,:),x2(r_idx,:));
    mask = abs(sum((x2_aug'.*(F*x1_aug')))) < threshold;
    num_in = sum(mask);
    if n_best < num_in
        n_best = num_in;
        y1 = x1(mask,:);
        y2 = x2(mask,:);
        idx = find(mask);
    end
end