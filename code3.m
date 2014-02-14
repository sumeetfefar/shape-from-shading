% Depth reconstruction

clear all


% Comment one of these to reconstruct Depth
% load('DataFile2.mat')
% load('DataFile1.mat')

load('DataFile3.mat')
load('DataFile4.mat')

M = size(E,1);
N = size(E,2);


zn = zeros(size(E));
z_old = zeros(size(E));

% Assume the Monge surface, calculate the first order partial derivatives 
% of the p and q
px = diff(pn,1,1);
qy = diff(qn,1,2);


% Apply the iterative method to estimate the value of depth at each point
iters = 1000;
for kk = 1:iters,
    disp(kk)
    for i=2:(M-1),
        for j=2:(N-1),
            if mask(i,j)==1
                zn(i,j) = Ravg(z_old,i,j) + px(i,j) + qy(i,j);
            else
                zn(i,j) = 0;
            end
        end
    end
    z_old = zn;
end

figure;
imshow(mat2gray(zn));