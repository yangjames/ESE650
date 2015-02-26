function params = baum_welch(O,M,N,max_iter)
T = length(O);

% initialize probabilities and transition matrix
Pi = zeros(N,1);
Pi(1) = 1;
A = diag(ones(1,N)*0.5) + circshift(diag(ones(1,N)*0.5),-1);
B = ones(N,M);
B = bsxfun(@rdivide, B, sum(B,2));

for iter = 1:max_iter

    % calculate alpha
    alpha_hat = zeros(N,T);
    Z_alpha = zeros(1,T);
    temp = B(:,1).*Pi;
    Z_alpha(1) = log(sum(temp));
    alpha_hat(:,1) = temp/exp(Z_alpha(1));
    for t = 2:T
        temp = (A'*alpha_hat(:,t-1)).*B(:,O(t));
        alpha_hat(:,t) = temp/max(sum(temp),eps);
        Z_alpha(t) = Z_alpha(t-1)+log(max(sum(temp),eps));
    end
    %{
    if sum(alpha_hat(:,end)*exp(Z_alpha(end))) > 0.99
        break;
    end
    %}
    % calculate beta
    beta_hat = zeros(N,T);
    Z_beta = zeros(1,T);
    temp = ones(N,1);
    Z_beta(end) = log(sum(temp));
    beta_hat(:,end) = temp/exp(Z_beta(end));
    for t = fliplr(1:T-1)
        temp = A*(B(:,O(t+1)).*beta_hat(:,t+1));
        beta_hat(:,t) = temp/max(sum(temp),eps);
        Z_beta(t) = Z_beta(t+1)+log(max(sum(temp),eps));
    end

    % calculate gamma
    %gamma = zeros(N,T);
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
end
params.A = A;
params.B = B;
params.Pi = Pi;
%params.C = C;