% Shape from Shading, EE 702, 2014
% Ashwin Kachhara, Sumeet Fefar

M=60;
N=60;

radius=25;

Depth = zeros(M,N);
E = 0.2 * ones(M,N);
% p = 0 * ones(M,N);
% q = 0 * ones(M,N);
% ptemp = 0 * ones(M,N);
% qtemp = 0 * ones(M,N);
% f = 0 * ones(M,N);
% g = 0 * ones(M,N);
s=[1,1];

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
        if(current_radius < radius)
            Depth(i,j) = round(sqrt(radius^2 - (i-M/2)^2 - (j-N/2)^2));
            p = (i-M/2)/Depth(i,j);
            q = (j-N/2)/Depth(i,j);
            mask(i,j)=1;

            temp = (s(1)*p + s(2)*q + 1)/sqrt( (s(2)^2 + s(1)^2 +1) * (p^2 + q^2 +1));
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
            
        end
    end
end


figure(1)
imshow(mat2gray(Depth));
imshow(mat2gray(E));

save('DataFile1.mat', 'E', 's','radius','mask','boundary');




