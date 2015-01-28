clear all
clc
%{d
%% load image data
fprintf('Loading data...\n')
load barrel_nobarrel_data.mat

%% convert barrel pixels to YCbCr and store
fprintf('Gathering pixel data...\n')
data = [];
for i = 1:length(file_names)
    im = rgb2ycbcr(imread(file_names{i}));
    Y = im(:,:,1);
    Cb = im(:,:,2);
    Cr = im(:,:,3);
    Y = Y(coords{1});
    Cb = Cb(coords{1});
    Cr = Cr(coords{1});

    clear im
    data = [data;Y Cb Cr];
    clear Y Cb Cr
end
data = double(data);
%}
%{
test_mu1 = [0.9 0.5 5];
test_sigma1 = 0.1*eye(3)
test_mu2 = [10 5 7];
test_sigma2 = 0.2*eye(3)
data = [mvnrnd(test_mu1,test_sigma1,1005); mvnrnd(test_mu2,test_sigma2,1005)];
%}
%{
figure(2)
clf
%plot(data(:,1),data(:,2),'.')
plot3(data(:,1),data(:,2),data(:,3),'.')
xlabel('Y')
ylabel('Cb')
zlabel('Cr')
grid on
axis equal
drawnow
%}
%% create gaussians and initialize GMM
fprintf('Initializing parameters...')
[n,d] = size(data);
k = 4;
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
done = 0;
gamma = zeros(n,k);
iterator = 1;
while ~done
    %fprintf('iteration: %d\n',iterator);
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
    
    % check exit condition
    if abs(L-L_new) < 0.001
        done = 1;
    end
    abs(L_new-L)
    L = L_new;
    iterator = iterator+1;
end
    
%% plot stuff
figure(1)
clf
plot3(data(1:100:end,1),data(1:100:end,2),data(1:100:end,3),'.')
axis equal
lims = [0 255];
grid on
hold on
for i = 1:k
    [U,L] = eig(A{i});
    N = 1;
    radii = N*sqrt(diag(L));
    [xc,yc,zc] = ellipsoid(0,0,0,radii(1),radii(2),radii(3));
    a = kron(U(:,1),xc); b = kron(U(:,2),yc); c = kron(U(:,3),zc);
    data = a+b+c; n = size(data,2);
    x = data(1:n,:)+mu(i,1); y = data(n+1:2*n,:)+mu(i,2); z = data(2*n+1:end,:)+mu(i,3);
    surf(x,y,z);
    alpha(0.1)
end
drawnow