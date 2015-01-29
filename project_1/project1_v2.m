clear all
clc
%{d
%% load image data
fprintf('Loading data...\n')
load barrel_nobarrel_data.mat

%% convert barrel pixels to YCbCr and store
fprintf('Gathering pixel data...\n')
data1 = [];
data2 = [];
for i = 1:length(file_names)
    % open image and convert to YCbCr
    im = rgb2ycbcr(imread(file_names{i}));
    
    % extract Y, Cb, and Cr channels
    Y = im(:,:,1);
    Cb = im(:,:,2);
    Cr = im(:,:,3);
    
    % get barrel colors
    Y1 = Y(coords{i});
    Cb1 = Cb(coords{i});
    Cr1 = Cr(coords{i});
    
    %{
    figure(1)
    clf
    mask = false(size(im,1),size(im,2));
    mask(coords{i}) = true;
    imshow(mask)
    pause
    %}
    
    % get non-barrel colors
    indices = (1:size(im,1)*size(im,2))';
    indices = indices(~ismember(indices,coords{i}));
    Y2 = Y(indices);
    Cb2 = Cb(indices);
    Cr2 = Cr(indices);
    
    %{
    figure(1)
    clf
    mask = true(size(im,1),size(im,2));
    mask(coords{i}) = false;
    imshow(mask)
    pause
    %}
    
    clear im
    data1 = [data1; Y1 Cb1 Cr1];
    data2 = [data2; Y2 Cb2 Cr2];
    clear Y Cb Cr Y1 Cb1 Cr1 Y2 Cb2 Cr2
end
data1 = double(data1);
data2 = double(data2);
%}
%{
%% intialize toy dataset
test_mu1 = [0.9 0.5 5];
test_sigma1 = 0.1*eye(3)
test_mu2 = [10 5 7];
test_sigma2 = 0.2*eye(3)
data = [mvnrnd(test_mu1,test_sigma1,1005); mvnrnd(test_mu2,test_sigma2,1005)];
%}
%{
figure(2)
clf
%plot(data1(:,1),data1(:,2),'.')
plot3(data1(:,1),data1(:,2),data1(:,3),'.')
xlabel('Y')
ylabel('Cb')
zlabel('Cr')
grid on
axis equal
drawnow
%}
%% create gaussians and initialize GMM
fprintf('Initializing parameters...\n')
[n,d] = size(data1);
k = 2;
A = cell(k,1);
mu = zeros(k,d);
w = zeros(1,k);
P = zeros(n,k);

for i = 1:k
    A{i} = diag(rand(d,1))*255;
    w(i) = 1/k;
    mu(i,:) = rand(1,d)*255;
    P(:,i) = compute_gaussian_density(data1,mu(i,:),A{i});
end
P(P<eps) = eps;

L = 1/n*sum(log(P*w'));

%% start EM
fprintf('Starting EM...\n')
done = 0;
gamma = zeros(n,k);
iterator = 1;
delta = 0.002;
L_new = 0;
current_delta = abs(L-L_new);
while ~done
    % progress print statement
    fprintf('iteration: %d | delta: %6.6f | target delta: %6.6f\n',iterator, current_delta, delta);
    
    % E-step
    gamma = bsxfun(@rdivide,bsxfun(@times,w,P),P*w');
    dumb_n = sum(gamma);
    
    % M-step
    w = dumb_n/n;
    for i = 1:k
        mu(i,:) = 1/dumb_n(i)*sum(bsxfun(@times,gamma(:,i),data1));
        mean_centered = bsxfun(@minus,data1,mu(i,:));
        A{i} = 1/dumb_n(i)*bsxfun(@times,gamma(:,i),mean_centered)'*mean_centered;
        P(:,i) = compute_gaussian_density(data1,mu(i,:),A{i});
    end
    P(P<eps) = eps;
    
    % calculate log likelihood
    L_new = 1/n*sum(log(P*w'));
    
    % check exit condition
    if abs(L-L_new) < delta
        done = 1;
    end
    current_delta = abs(L-L_new);
    L = L_new;
    iterator = iterator+1;
end
model1.weight = w;
model1.mean = mu;
model1.cov = A;
model1.num_clusters = k;

%{d
%% plot stuff
figure(1)
clf
sparsity = 100;
plot(data1(1:sparsity:end,1),data1(1:sparsity:end,2),'.')
axis equal
grid on
hold on
for i = 1:model1.num_clusters
    
end
%}
%{d
%% plot stuff
figure(1)
clf
sparsity = 100;
plot3(data1(1:sparsity:end,1),data1(1:sparsity:end,2),data1(1:sparsity:end,3),'.')
axis equal
lims = [0 255];
grid on
hold on
xlabel('Y')
ylabel('Cb')
zlabel('Cr')
for i = 1:k
    % get principle components of covariance matrix
    [V,Lambda] = eig(model1.cov{i});
    R = 1;
    prin_comps = R*sqrt(diag(Lambda));
    
    % make an ellipsoid
    [xe,ye,ze] = ellipsoid(0,0,0,prin_comps(1),prin_comps(2),prin_comps(3));
    
    % rotate axis aligned ellipse
    rot_el = kron(V(:,1),xe) + kron(V(:,2),ye) + kron(V(:,3),ze);
    
    % get ellipse surface points
    rot_x = rot_el(1:size(rot_el,2),:)+model1.mean(i,1);
    rot_y = rot_el(size(rot_el,2)+1:2*size(rot_el,2),:)+model1.mean(i,2);
    rot_z = rot_el(2*size(rot_el,2)+1:end,:)+model1.mean(i,3);
    
    % plot ellipse
    surf(rot_x,rot_y,rot_z);
    alpha(0.1)
end
drawnow
%}

%% test classifications
for i = 1:length(file_names)
    % get image
    rgbim = imread(file_names{i});
    G = fspecial('gaussian',[5 5],4);
    rgbim = imfilter(rgbim,G,'same');
    
    im = rgb2ycbcr(rgbim);
    [r,c,~] = size(im);
    % get pixel values
    Y = reshape(im(:,:,1),r*c,1);
    Cb = reshape(im(:,:,2),r*c,1);
    Cr = reshape(im(:,:,3),r*c,1);
    colors = double([Y Cb Cr]);
    
    % calculate probability density and Mahalanobis distance of each pixel
    P = zeros(r*c,model1.num_clusters);
    DM = zeros(r*c,model1.num_clusters);
    for j = 1:model1.num_clusters
        P(:,j) = compute_gaussian_density(colors,model1.mean(j,:),model1.cov{j});
        mean_centered = bsxfun(@minus,colors,model1.mean(j,:));
        DM(:,j) = sqrt(sum((mean_centered/A{j}).*mean_centered,2));
    end
    
    % check which pixels fall close to a cluster
    cdm = DM>6;
    
    % create mask for pixels
    pre_filt_mask = reshape(uint8(~(sum(cdm,2) > 0))*255,r,c);
    mask = medfilt2(pre_filt_mask);
        
    % display image
    figure(2)
    clf
    imshow(rgbim)
    drawnow
    
    % display mask
    figure(3)
    clf
    imshow(mask)
    drawnow
    
    % wait until a key is pressed
    pause
end