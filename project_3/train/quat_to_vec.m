function v = quat_to_vec(q)

alpha = 2*acos(q(1,:));
e_hat = q(2:end,:);
v = bsxfun(@times,alpha./sin(alpha/2),e_hat);
v(:,alpha == 0) = repmat([0 0 0]',1,sum(alpha==0));