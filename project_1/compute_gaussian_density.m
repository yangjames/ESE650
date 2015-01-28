function P = compute_gaussian_density(data,mu,A)
centered_data = bsxfun(@minus,data,mu);
P = 1/( (2*pi)^(size(data,2)/2) * sqrt(det(A)))...
    *exp( -1/2 * sum( (centered_data/A) .* centered_data ,2));
end