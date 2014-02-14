% Mozart bust depth map provided by R. Zhang, P.-S Tsai, J.E. Cryer 
% and M. Shah. (Computer vision Lab. of UCF)

% load('moz256.mat')
% moz256 = imresize(moz256, 0.25);

moz256 = rgb2gray(imread('mozart.jpg'));


[M,N] = size(moz256);

E = moz256;

s_orig=[0,0];

p_init = zeros(M,N);
q_init = zeros(M,N);
f_init = zeros(M,N);
g_init = zeros(M,N);

mask = zeros(M,N);

boundary = zeros(M,N);



for i=2:M-1,
    for j=2:N-1,
       if(Ravg(moz256,i,j)>5)
           mask(i,j)=1;
       end
            

    end
end

mask = medfilt2(mask);
ed = edge(mask);

mask = mask | ed ;
% figure
% imshow(mat2gray(mask));
% figure
% imshow(mat2gray(ed));

% [Gmag, Gdir] = imgradient(mask);
% 
% f = 2*cos(Gdir) ;
% g = 2*sin(Gdir) ;

for i=2:M-1,
    for j=2:N-1,
        ftemp = (2*mask(i+1,j) + mask(i+1,j+1) + mask(i+1,j-1)) - (2*mask(i-1,j) + mask(i-1,j+1) +mask(i-1,j-1));
        gtemp = (2*mask(i,j+1) + mask(i+1,j+1) + mask(i-1,j+1)) - (2*mask(i,j-1) + mask(i-1,j-1) +mask(i+1,j-1));
        if(gtemp ~=0 || ftemp ~=0)
            f_init(i,j)= - 2 *ftemp * mask(i,j)/sqrt(ftemp^2+gtemp^2);
            g_init(i,j)= - 2 *gtemp * mask(i,j)/sqrt(ftemp^2+gtemp^2);
            boundary(i,j)=1;
        end
    end
end

% figure
% imshow(mat2gray(p_init));
% figure
% imshow(mat2gray(q_init));


%Source passed to algorithm used for reconstruction
s = [0,1];
% 
figure;
imshow(mat2gray(E));


En = double(E)/255;

p_old = f_init;
q_old = g_init;
pn = zeros(size(En));
qn = zeros(size(En));

iters = 500;
lambda = 1;



for kk = 1:iters,
    disp(kk)
    for i=2:(M-1),
        for j=2:(N-1),
            if(boundary(i,j)==0 && mask(i,j) ==1)
                pn(i,j) = Ravg(p_old,i,j) + (1/lambda)*( En(i,j) - Rfgval(p_old(i,j), q_old(i,j), s))*Rf(p_old(i,j),q_old(i,j),s(1),s(2));
                qn(i,j) = Ravg(q_old,i,j) + (1/lambda)*( En(i,j) - Rfgval(p_old(i,j), q_old(i,j), s))*Rg(p_old(i,j),q_old(i,j),s(1),s(2));
            else 
                pn(i,j) = p_old(i,j);
                qn(i,j) = q_old(i,j);
            end
        end
    end
    for i=2:M-1
        pn(i,N)=pn(i,N-1);
        pn(i,1)=pn(i,2);
        qn(i,N)=qn(i,N-1);
        qn(i,1)=qn(i,2);
    end
    for i=2:N-1
        pn(N,i)=pn(N-1,i);
        pn(1,i)=pn(2,i);
        qn(N,i)=qn(N-1,i);
        qn(1,i)=qn(2,i);
    end
    
    pn(1,N)=pn(2,N-1);   
    pn(N,1)=pn(N-1,2);   
    pn(1,1)=pn(2,2);   
    pn(N,N)=pn(N-1,N-1);   
    
    qn(1,N)=qn(2,N-1);   
    qn(N,1)=qn(N-1,2);   
    qn(1,1)=qn(2,2);   
    qn(N,N)=qn(N-1,N-1);   
    
    p_old = pn;
    q_old = qn;
end

% figure;
% imshow(mat2gray(qn));


zn = zeros(size(E));
z_old = zeros(size(E));

px = diff(pn,1,1);
qy = diff(qn,1,2);


iters2=1000;

for kk = 1:iters2,
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
    for i=2:M-1
        zn(i,N)=zn(i,N-1);
        zn(i,1)=zn(i,2);
    end
    for i=2:N-1
        zn(N,i)=zn(N-1,i);
        zn(1,i)=zn(2,i);
    end
    
    zn(1,N)=zn(2,N-1);   
    zn(N,1)=zn(N-1,2);   
    zn(1,1)=zn(2,2);   
    zn(N,N)=zn(N-1,N-1);   
    
    z_old = zn;
end

figure;
surf(mat2gray(zn));

