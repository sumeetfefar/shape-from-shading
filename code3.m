clear all

load('DataFile1.mat')
load('DataFile2.mat')

M = size(R_est,1);
N = size(R_est,2);

iters = 1000;

zn = zeros(size(R_est));
zold = zeros(size(R_est));

px = diff(pn,1,1);
qy = diff(qn,1,2);

for kk = 1:iters,
    disp(kk)
    for i=2:(M-1),
        for j=2:(N-1),
            if mask(i,j)==1
                zn(i,j) = (zold(i,j-1)+zold(i-1,j)+zold(i,j+1)+zold(i+1,j))/4 + px(i,j) + qy(i,j);
            else
                zn(i,j) = 0;
            end
        end
    end
    zold = zn;
end

figure;
imshow(mat2gray(zn));