% Double Source with weight a & b code Reconstruction


load('DataFile3.mat');

M = size(E,1);
N = size(E,2);


% En = E_noise;
En = E;

p_old = p_init;
q_old = q_init;
pn = zeros(size(En));
qn = zeros(size(En));




iters = 300;
lambda = 0.5  ;

for kk = 1:iters,
    disp(kk)
    for i=2:(M-1),
        for j=2:(N-1),
            if(boundary(i,j)==0 && mask(i,j) ==1)
                pn(i,j) = Ravg(p_old,i,j) + (1/lambda)*( En(i,j) - (a*Rval(p_old(i,j),q_old(i,j),s1)+b*Rval(p_old(i,j),q_old(i,j),s2)))*(a*Rp(p_old(i,j),q_old(i,j),s1(1),s1(2))+b*Rp(p_old(i,j),q_old(i,j),s2(1),s2(2)));
                qn(i,j) = Ravg(q_old,i,j) + (1/lambda)*( En(i,j) - (a*Rval(p_old(i,j),q_old(i,j),s1)+b*Rval(p_old(i,j),q_old(i,j),s2)))*(a*Rq(p_old(i,j),q_old(i,j),s1(1),s1(2))+b*Rq(p_old(i,j),q_old(i,j),s2(1),s2(2)));
            else 
                pn(i,j) = p_old(i,j);
                qn(i,j) = q_old(i,j);
            end
        end
    end
    p_old = pn;
    q_old = qn;
end

save('DataFile3.mat', 'pn', 'qn');