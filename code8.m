% Mozart bust depth map provided by R. Zhang, P.-S Tsai, J.E. Cryer 
% and M. Shah. (Computer vision Lab. of UCF)

load('moz256.mat')
moz256 = imresize(moz256, 0.25);

radius = 10;

[M,N] = size(moz256);

E = 0.2 * ones(M,N);

s_orig=[2,1];

p_init = zeros(M,N);
q_init = zeros(M,N);
f_init = zeros(M,N);
g_init = zeros(M,N);
p_orig = zeros(M,N);
q_orig = zeros(M,N);
f_orig = zeros(M,N);
g_orig = zeros(M,N);
mask = zeros(M,N);
boundary = zeros(M,N);

for i=1:M,
    for j=1:N,
        current_radius = sqrt((i-M/2)^2 + (j-N/2)^2);
%         if(current_radius < radius)
            p = (i-M/2)/moz256(i,j);
            q = (j-N/2)/moz256(i,j);
            mask(i,j)=1;

            temp = Rval(p, q, s_orig);
            if (round(current_radius) == round(radius))
                p_init(i,j)= p;
                q_init(i,j)= q;
                boundary(i,j)=1;
            end
            p_orig(i,j) = p;
            q_orig(i,j) = q;
            if (round(current_radius) >= radius -1 && round(current_radius) < radius+1)
                f_init(i,j) = 2*p/(1+sqrt(1+p^2+q^2));
                g_init(i,j) = 2*q/(1+sqrt(1+p^2+q^2));
            end
            f_orig(i,j) = 2*p/(1+sqrt(1+p^2+q^2));
            g_orig(i,j) = 2*q/(1+sqrt(1+p^2+q^2));
            if(temp>0)
                E(i,j) = temp;
            else
                E(i,j) = 0;

            end
            
%         end
    end
end

%Source passed to algorithm used for reconstruction
s = [0,0];
% 
% figure;
% imshow(mat2gray(E));

% Noisy image with gaussian noise
E_noise = imnoise(E,'gaussian',0,5);

En = E;

p_old = p_init;
q_old = q_init;
pn = zeros(size(En));
qn = zeros(size(En));

iters = 700;
lambda = 0.02;

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

zn = zeros(size(E));
z_old = zeros(size(E));

px = diff(pn,1,1);
qy = diff(qn,1,2);

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
surf(mat2gray(zn));

