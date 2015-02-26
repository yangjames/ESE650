function params = baum_welch(O,N,max_iter)
T = length(O);
M = length(unique(O));

% initialize probabilities and transition matrix
Pi = zeros(N,1);
Pi(1) = 1;
%Pi = ones(N,1)/N;
A = diag(ones(1,N)*0.5) + circshift(diag(ones(1,N)*0.5),-1);
B = ones(N,M);
B = bsxfun(@rdivide, B, sum(B,2));

A_prev = A;
B_prev = B;
Pi_prev = Pi;

logP = inf(1,max_iter);
for iter = 1:max_iter
    
    % calculate alpha
    alpha_hat = zeros(N,T);
    Z_alpha = zeros(1,T);
    temp = B(:,1).*Pi;
    Z_alpha(1) = log(sum(temp));
    alpha_hat(:,1) = temp/max(exp(Z_alpha(1)),eps);
    for t = 2:T
        temp = (A'*alpha_hat(:,t-1)).*B(:,O(t));
        alpha_hat(:,t) = temp/max(sum(temp),eps);
        Z_alpha(t) = Z_alpha(t-1)+log(max(sum(temp),eps));
    end
    
    % calculate beta
    beta_hat = zeros(N,T);
    Z_beta = zeros(1,T);
    temp = ones(N,1);
    Z_beta(end) = log(sum(temp));
    beta_hat(:,end) = temp/max(exp(Z_beta(end)),eps);
    for t = fliplr(1:T-1)
        temp = A*(B(:,O(t+1)).*beta_hat(:,t+1));
        beta_hat(:,t) = temp/max(sum(temp),eps);
        Z_beta(t) = Z_beta(t+1)+log(max(sum(temp),eps));
    end

    % calculate gamma
    gamma = bsxfun(@rdivide,alpha_hat.*beta_hat,max(sum(alpha_hat.*beta_hat),eps));

    % calculate ksi
    ksi = zeros(N,N,T-1);
    ksi_hat = zeros(N,N,T-1);
    Z_ksi = zeros(1,T-1);
    for t = 1:T-1
        ksi_hat(:,:,t) = bsxfun(@times,bsxfun(@times,alpha_hat(:,t),A),(B(:,O(t+1)).*beta_hat(:,t+1))');
        Z_ksi(t) = log(sum(sum(ksi_hat(:,:,t))));
        ksi(:,:,t) = ksi_hat(:,:,t)/max(exp(Z_ksi(t)),eps);
    end

    % calculate new model parameters
    Pi = gamma(:,1);
    A = bsxfun(@rdivide,sum(ksi,3),max(sum(gamma(:,1:T-1),2),eps));
    obs = unique(O);
    B = zeros(size(B));
    for i = 1:length(obs)
        B(:,obs(i)) = sum(gamma(:,O==obs(i)),2)./max(sum(gamma,2),eps);
    end
    B(B<eps) = eps;
    B = bsxfun(@rdivide,B,sum(B,2));
    
    % Viterbi algorithm initialization
    phi = zeros(N,T);
    phi(:,1) = log(Pi)+log(B(:,O(1)));
    
    % calculate probabilities
    for t = 2:T
        phi(:,t) = max(bsxfun(@plus,phi(:,t-1),log(A)))'...
                + log(B(:,O(t)));
    end
    logP(iter) = max(phi(:,end));
    
    if iter > 1
        if abs(logP(iter)-logP(iter-1)) < 0.005
            A = A_prev;
            B = B_prev;
            Pi = Pi_prev;
            fprintf('lol\n')
            break
        else
            A_prev = A;
            B_prev = B;
            Pi_prev = Pi;
        end
    end
    
    % find the winner
    if iter > 1
        fprintf('delta logprob: %6.6f | logprob: %6.6f\n', logP(iter)-logP(iter-1), logP(iter))
    end
end
params.A = A;
params.B = B;
params.Pi = Pi;