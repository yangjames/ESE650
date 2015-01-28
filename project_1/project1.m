clear all
clc

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

%% create gaussians and initialize GMM
fprintf('Initializing parameters...\n')
k = 4;
A = cell(k,1);
mu = zeros(k,3);
w = zeros(1,k);
for i = 1:k
    A{i} = eye(3);
    w(i) = 1/k;
    mu(i,:) = rand(1,3)*255;
end

%% start EM
P = zeros(size(data,1),k);
max_iteration = 200;
done = 0;
alpha = zeros(1,k);

% initialize gaussians
fprintf('Initializing gaussians...\n')
for i = 1:k
    P(:,i) = compute_gaussian_density(data,mu(i,:),A{i});
end
P(P<eps) = eps;

% calculate initial log likelihood
L = 1/size(data,1)*sum(log(sum(P,2)))

% continue EM algorithm
fprintf('Staring EM...\n')
iteration = 1;
while ~done %iteration < max_iteration || done
    fprintf('iteration: %d\n',iteration)
    
    % E-step
    alpha = sum(P)/size(data,1);
    
    % M-step
    for i = 1:k
        mu(i,:) = mean(alpha(i)*data);
        mean_centered = bsxfun(@minus, data,mu(i,:));
        A{i} = 1/size(data,1)*((alpha(i)*mean_centered)'*mean_centered);
    end
   
    % compute gaussians
    for i = 1:k
        P(:,i) = compute_gaussian_density(data,mu(i,:),A{i});
    end
    P(P<eps) = eps;
    
    % calculate log likelihood
    L_new = 1/size(data,1)*sum(log(sum(P,2)));
    
    % check exit condition
    if abs(L_new - L) > 0.1
        done = 1;
    end
    L_new - L
    L = L_new;
    iteration = iteration+1;
end

%% plot stuff
figure(1)
clf
plot3(data(1:100:end,1),data(1:100:end,2),data(1:100:end,3),'.')
axis equal
lims = [0 255];
xlim(lims)
ylim(lims)
zlim(lims)
grid on
hold on
for i = 1:k
    [U,L] = eig(A{i});
    N = 1;
    radii = N*sqrt(diag(L));
    [xc,yc,zc] = ellipsoid(0,0,0,radii(1),radii(2),radii(3));
    a = kron(U(:,1),xc); b = kron(U(:,2),yc); c = kron(U(:,3),zc);
    data = a+b+c; n = size(data,2);
    x = data(1:n,:)+M(1); y = data(n+1:2*n,:)+M(2); z = data(2*n+1:end,:)+M(3);
    sc = surf(x,y,z);
    %axis equal
    alpha(0.5)
end
