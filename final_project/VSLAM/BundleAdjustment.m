function [Cr_set, Rr_set, X3D] = BundleAdjustment(K, Cr_set, Rr_set, X3D, ReconX, Mx_bundle, My_bundle)

cP = cell(size(Rr_set));

mask = logical(ReconX);
num_features = sum(mask);
num_frames = length(Rr_set);

X = X3D(mask,:);
Mx = Mx_bundle(mask,:);
My = My_bundle(mask,:);
%V = V_bundle(mask,:);
measurements = zeros(2*num_features,num_frames);
for i = 1:num_frames
    cP{i} = K*Rr_set{i}*[eye(3) -Cr_set{i}];
    measurements(:,i) = reshape([Mx(:,i) My(:,i)]',num_features*2,1);
end

[cP, X] = sba_wrapper(measurements, cP, X, K);
X3D(mask,:) = X;
for i = 1:num_frames
    H = K\cP{i};
    Rr_set{i} = H(1:3,1:3);
    Cr_set{i} = H(:,4);
end