% p,q estimate from single source

%clear all

% load('DataFile1.mat');
% load('PQFfile.mat');

M = size(E,1);
N = size(E,2);


% En = E_noise;
En = E;

p_old = p_init;
q_old = q_init;
pn = zeros(size(En));
qn = zeros(size(En));




iters = 300;
lambda = 0.75  ;

for kk = 1:iters,
    disp(kk)
    for i=2:(M-1),
        for j=2:(N-1),
            if(boundary(i,j)==0 && mask(i,j) ==1)
                pn(i,j) = Ravg(p_old,i,j) + (1/lambda)*( En(i,j) - Rval(p_old(i,j), q_old(i,j), s))*Rp(p_old(i,j),q_old(i,j),s(1),s(2));
                qn(i,j) = Ravg(q_old,i,j) + (1/lambda)*( En(i,j) - Rval(p_old(i,j), q_old(i,j), s))*Rq(p_old(i,j),q_old(i,j),s(1),s(2));
            else 
                pn(i,j) = p_old(i,j);
                qn(i,j) = q_old(i,j);
            end
        end
    end
    p_old = pn;
    q_old = qn;
end

save('DataFile2.mat', 'pn', 'qn');
