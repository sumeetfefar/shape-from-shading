clear all

load('DataFile1.mat');
% load('PQFfile.mat');

p_old= p_init;
q_old = q_init;
pn = zeros(size(E));
qn = zeros(size(E));


M = size(E,1);
N = size(E,2);

iters = 1000;
lambda = 50  ;

for kk = 1:iters,
    disp(kk)
    for i=2:(M-1),
        for j=2:(N-1),
            if(boundary(i,j)==0 && mask(i,j) ==1)
                pn(i,j) = (p_old(i,j-1)+p_old(i-1,j)+p_old(i,j+1)+p_old(i+1,j))/4 + (1/lambda)*(E(i,j) - Rval(p_old(i,j), q_old(i,j), s))*Rp(p_old(i,j),q_old(i,j),s(1),s(2));
                qn(i,j) = (q_old(i,j-1)+q_old(i-1,j)+q_old(i,j+1)+q_old(i+1,j))/4 + (1/lambda)*(E(i,j) - Rval(p_old(i,j), q_old(i,j), s))*Rq(p_old(i,j),q_old(i,j),s(1),s(2));
            else
                pn(i,j) = p_old(i,j);
                qn(i,j) = q_old(i,j);
            end
        end
    end
    p_old = pn;
    q_old = qn;
end

R_est = zeros(size(E));
for i=1:M,
    for j=1:N,
        if(mask(i,j)==1)
            R_est(i,j) = Rval(pn(i,j), qn(i,j), s);
        end
    end
end
figure(2)
imshow(mat2gray(R_est));

save('DataFile2.mat', 'R_est', 'pn', 'qn');
