function q_bar = quat_conj(q)
q_bar = bsxfun(@times,q,[1;-1;-1;-1]);