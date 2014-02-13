clear all

load('DataFile1.mat');

p_old= p_init;
q_old = q_init;
pn = zeros(size(E));
qn = zeros(size(E));

M = size(E,1);
N = size(E,2);

iters = 1000;
lambda = 0;
mse = zeros(size((105:5:150),2),1);

for lambda = 105:5:150,
    disp(lambda/5)
    for kk = 1:iters,
%         disp(kk)
        for i=2:(M-1),
            for j=2:(N-1),
                if(boundary(i,j)==0 && mask(i,j) ==1)
                    pn(i,j) = Ravg(p_old,i,j) + (1/lambda)*(E(i,j) - Rval(p_old(i,j), q_old(i,j), s))*Rp(p_old(i,j),q_old(i,j),s(1),s(2));
                    qn(i,j) = Ravg(q_old,i,j) + (1/lambda)*(E(i,j) - Rval(p_old(i,j), q_old(i,j), s))*Rq(p_old(i,j),q_old(i,j),s(1),s(2));
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
%     figure(2)
%     imshow(mat2gray(R_est));
    imwrite(mat2gray(pn), strcat('img/pn-',num2str(lambda/5),'.jpg'));
    imwrite(mat2gray(qn), strcat('img/qn-',num2str(lambda/5),'.jpg'));

    zn = zeros(size(R_est));
    zold = zeros(size(R_est));

    px = diff(pn,1,1);
    qy = diff(qn,1,2);

    for kk = 1:iters,
%         disp(kk)
        for i=2:(M-1),
            for j=2:(N-1),
                if mask(i,j)==1
                    zn(i,j) = Ravg(zold,i,j) + px(i,j) + qy(i,j);
                else
                    zn(i,j) = 0;
                end
            end
        end
        zold = zn;
    end

%     figure(3);
%     imshow(mat2gray(zn));
    imwrite(mat2gray(zn), strcat('img/zn-',num2str(lambda/5),'.jpg'));

    mse(lambda/5) = sum(sum((zn-Depth).^2))/(M*N);
end
