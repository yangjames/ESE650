function [Cset, Rset] = ExtractCameraPose(E)
Cset = cell(4,1);
Rset = cell(4,1);

[U,~,V] = svd(E);
W = [0 -1 0; 1 0 0; 0 0 1];

Rset{1} = U*W*V';
Cset{1} = U(:,3);
if det(Rset{1}) < 0
    Rset{1} = -Rset{1};
    Cset{1} = -Cset{1};
end

Rset{2} = U*W*V';
Cset{2} = -U(:,3);
if det(Rset{2}) < 0
    Rset{2} = -Rset{2};
    Cset{2} = -Cset{2};
end

Rset{3} = U*W'*V';
Cset{3} = U(:,3);
if det(Rset{3}) < 0
    Rset{3} = -Rset{3};
    Cset{3} = -Cset{3};
end

Rset{4} = U*W'*V';
Cset{4} = -U(:,3);
if det(Rset{4}) < 0
    Rset{4} = -Rset{4};
    Cset{4} = -Cset{4};
end