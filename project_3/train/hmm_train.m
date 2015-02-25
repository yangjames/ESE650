clear all
close all
clc

pattern = {'beat3','beat4','circle','eight','inf','wave'};
contents = dir;
file_names = cell(length(contents),1);
for file_idx = 1:length(contents)
    file_names{file_idx} = contents(file_idx).name;...
end

for pattern_num = 1:length(pattern)
    acc = [];
    orient = [];
    time = [];
    valid_files = find(~cellfun(@isempty,regexp(file_names,['^' pattern{pattern_num} '.+(\.mat)$'])));
    for i = 1:length(valid_files)
        data = load(file_names{valid_files(i)});
        acc = [acc data.full_data.acceleration/9.81];
        orient = [orient data.full_data.orientation];
        time = [time data.full_data.time'];
    end
    
    acc_g = zeros(size(acc));
    for i = 1:length(time)
        acc_i = quat_mult(quat_mult(orient(:,i),[0; acc(:,i)]),quat_conj(orient(:,i)));
        acc_g(:,i) = acc_i(2:4);
    end
    
    N = 10;
    M = 8;
    [idx,C] = kmeans(orient',M);
    %idx = repmat([1 2 3 4],1,20);
    T = length(idx);
    
    % initialize probabilities and transition matrix
    Pi = zeros(N,1);%ones(N,1)/N;
    Pi(1) = 1;
    A = diag(ones(1,N)*0.5) + circshift(diag(ones(1,N)*0.5),-1);
    B = rand(N,M);%ones(N,M);
    B = bsxfun(@rdivide, B, sum(B,2));
    
    % calculate alpha
    alpha = zeros(N,T);
    alpha(:,1) = B(:,1).*Pi;
    
    alpha_hat = zeros(N,T);
    Z_alpha = zeros(1,T);
    temp = B(:,1).*Pi;
    Z_alpha(1) = log(sum(temp));
    alpha_hat(:,1) = temp/exp(Z_alpha(1));
    for t = 2:T
        alpha(:,t) = (A'*alpha(:,t-1)).*B(:,idx(t));
        
        temp = (A'*alpha_hat(:,t-1)).*B(:,idx(t));
        alpha_hat(:,t) = temp/sum(temp);
        Z_alpha(t) = Z_alpha(t-1)+log(sum(temp));
    end
    
    % calculate beta
    beta = zeros(N,T);
    beta(:,end) = 1;
    
    beta_hat = zeros(N,T);
    Z_beta = zeros(1,T);
    temp = ones(N,1);
    Z_beta(end) = log(sum(temp));
    beta_hat(:,end) = temp/exp(Z_beta(end));
    for t = fliplr(1:T-1)
        beta(:,t) = (A*B(:,idx(t+1))).*beta(:,t+1);
        
        %temp = (A*B(:,idx(t+1))).*beta_hat(:,t+1);
        temp = A*(B(:,idx(t+1)).*beta_hat(:,t+1));
        beta_hat(:,t) = temp/sum(temp);
        Z_beta(t) = Z_beta(t+1)+log(sum(temp));
    end
    
    % calculate gamma
    gamma = zeros(N,T);
    gamma = bsxfun(@rdivide,alpha_hat.*beta_hat,sum(alpha_hat.*beta_hat));
    
    % calculate ksi
    ksi = zeros(N,N,T-1);
    ksi_hat = zeros(N,N,T-1);
    Z_ksi = zeros(1,T-1);
    for t = 1:T-1
        ksi_hat(:,:,t) = bsxfun(@times,bsxfun(@times,alpha_hat(:,t),A'),(B(:,idx(t+1)).*beta_hat(:,t+1))');
        Z_ksi(t) = log(sum(sum(ksi_hat(:,:,t))));
        ksi(:,:,t) = ksi_hat(:,:,t)/exp(Z_ksi(t));
    end
    
    % calculate new model parameters
    Pi_bar = gamma(:,1);
    A_bar = bsxfun(@rdivide,sum(ksi,3),sum(gamma(:,1:T-1),2));
    obs = unique(idx);
    B_bar = zeros(size(B));
    for i = 1:length(obs)
        B_bar(:,obs(i)) = sum(gamma(:,idx==obs(i)),2)./sum(gamma,2);
    end
    params.A = A_bar;
    params.B = B_bar;
    params.Pi = Pi_bar;
    save([pattern{pattern_num} '_params.mat'],'params')
end