% p,q estimate from single source

clear all

load('DataFile1.mat');
% load('PQFfile.mat');

M = size(E,1);
N = size(E,2);

% s=[ double(s(1)) double(s(2))];
s=[ double(2*s(1)/(1+sqrt(1+s(1)^2+s(2)^2))) double(2*s(2)/(1+sqrt(1+s(1)^2+s(2)^2)))];

En = E_noise;
% En = E;

f_old = f_init;
g_old = g_init;
fn = zeros(size(En));
gn = zeros(size(En));

iters = 500;
lambda = 01  ;

for kk = 1:iters,
%     disp(kk)
    for i=2:(M-1),
        for j=2:(N-1),
            if(boundaryfg(i,j)==0 && mask(i,j) ==1)
                fn(i,j) = Ravg(f_old,i,j) + (1/lambda) * ( En(i,j) - Rfgval(f_old(i,j),g_old(i,j),s) ) * Rf(f_old(i,j),g_old(i,j),s(1),s(2));
                gn(i,j) = Ravg(g_old,i,j) + (1/lambda) * ( En(i,j) - Rfgval(f_old(i,j),g_old(i,j),s) ) * Rg(f_old(i,j),g_old(i,j),s(1),s(2));
            else 
                fn(i,j) = f_old(i,j);
                gn(i,j) = g_old(i,j);
            end
        end
    end
    f_old = fn;
    g_old = gn;
end

zn = zeros(size(E));
z_old = zeros(size(E));

fx = diff(fn,1,1);
fy = diff(fn,1,2);
gx = diff(gn,1,1);
gy = diff(gn,1,2);

for kk = 1:iters,
%     disp(kk)
    for i=2:(M-1),
        for j=2:(N-1),
            if mask(i,j)==1
                zn(i,j) = Ravg(z_old,i,j) + fx(i,j) + gy(i,j);% 4*((fx(i,j)+gy(i,j))*(4-fn(i,j)^2-gn(i,j)^2)+2*fn(i,j)*(fn(i,j)*fx(i,j)+gn(i,j)*gx(i,j))+2*gn(i,j)*(fn(i,j)*fy(i,j)+gn(i,j)*gy(i,j)))/(4-fn(i,j)^2-gn(i,j)^2)^2;
            else
                zn(i,j) = 0;
            end
        end
    end
    z_old = zn;
end

figure;
imshow(mat2gray(zn));
