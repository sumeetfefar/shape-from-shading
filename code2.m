load('Esfile.mat');
% load('PQFfile.mat');

p = zeros(size(E));
q = zeros(size(E));
ptemp = zeros(size(E));
qtemp = zeros(size(E));

M = size(E,1);
N = size(E,2);

iters = 10;
lamdba = 0.5;

for kk = 1:iters,
    ptemp = p;
    qtemp = q;
    for i=1:M,
        for j=1:N,
            p(i,j) = (ptemp(i,j-1)+ptemp(i-1,j)+ptemp(i,j+1)+ptemp(i+1,j))/4 + (1/lambda)*(E(i,j) - Rval(ptemp(i,j), qtemp(i,j), s))*(Rval(ptemp(i+1,j), qtemp(i,j), s)-Rval(ptemp(i-1,j), qtemp(i,j), s))/2;
            q(i,j) = (qtemp(i,j-1)+qtemp(i-1,j)+qtemp(i,j+1)+qtemp(i+1,j))/4 + (1/lambda)*(E(i,j) - Rval(ptemp(i,j), qtemp(i,j), s))*(Rval(ptemp(i,j), qtemp(i,j+1), s)-Rval(ptemp(i,j), qtemp(i,j-1), s))/2;
        end
    end
end

R_est = zeros(size(E));
for i=1:M,
    for j=1:N,
        R_est(i,j) = Rval(p(i,j), q(i,j), s);
    end
end

imshow(mat2gray(R_est));