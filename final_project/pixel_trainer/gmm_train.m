function model = gmm_train(data,k, delta)

%% create gaussians and initialize GMM
fprintf('Initializing parameters...\n')
[n,d] = size(data);
A = cell(k,1);
mu = zeros(k,d);
w = zeros(1,k);
P = zeros(n,k);

for i = 1:k
    A{i} = diag(rand(d,1))*255;
    w(i) = 1/k;
    mu(i,:) = rand(1,d)*255;
    P(:,i) = compute_gaussian_density(data,mu(i,:),A{i});
end
P(P<eps) = eps;

L = 1/n*sum(log(P*w'));

%% start EM
fprintf('Starting EM...\n')
done = 0;
iterator = 1;
L_new = 0;
while ~done
    
    % E-step
    gamma = bsxfun(@rdivide,bsxfun(@times,w,P),P*w');
    dumb_n = sum(gamma);
    
    % M-step
    w = dumb_n/n;
    for i = 1:k
        mu(i,:) = 1/dumb_n(i)*sum(bsxfun(@times,gamma(:,i),data));
        mean_centered = bsxfun(@minus,data,mu(i,:));
        A{i} = 1/dumb_n(i)*bsxfun(@times,gamma(:,i),mean_centered)'*mean_centered;
        P(:,i) = compute_gaussian_density(data,mu(i,:),A{i});
    end
    P(P<eps) = eps;
    
    % calculate log likelihood
    L_new = 1/n*sum(log(P*w'));
    current_delta = abs(L-L_new);
    
    % progress print statement
    fprintf('iteration: %d | delta: %6.6f | target delta: %6.6f\n',iterator, current_delta, delta);
    
    % check exit condition
    if current_delta < delta
        done = 1;
    end
    
    
    L = L_new;
    iterator = iterator+1;
end
model.weight = w;
model.mean = mu;
model.cov = A;
model.num_clusters = k;