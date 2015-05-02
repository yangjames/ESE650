function [C, R, X] = DisambiguateCameraPose(Cset,Rset,Xset)

num_poses = length(Cset);
counts = zeros(4,1);
for i = 1:num_poses
    counts(i) = sum(Rset{i}(3,:)*bsxfun(@minus,Xset{i}',Cset{i}) > 0);
end
[~,idx] = max(counts);
C = Cset{idx};
R = Rset{idx};
X = Xset{idx};