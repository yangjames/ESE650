load depth0.mat
load rgb0.mat
addpath('mex')
DEPTH_MAX = 4500; 
DEPTH_MIN = 400;

for k=2:2:numel(DEPTH)
   
   D = DEPTH{k}.depth';
   D = flip(D,2);
   D(D(:) <= DEPTH_MIN) = 0;
   D(D(:) >= DEPTH_MAX) = 0;  
   figure(1), imagesc(D); 
   
   
   R = djpeg(RGB{k}.jpeg);
   R = flip(R,2);
   figure(2), imagesc(R); 
   
   pause(0.25);
   
end

