% Shape from Shading, EE 702, 2014
% Ashwin Kachhara, Sumeet Fefar

% IMage generator with source direction = s_orig

M=64;
N=64;

radius=25;

Depth = zeros(M,N);
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
        if(current_radius < radius)
            Depth(i,j) = round(sqrt(radius^2 - (i-M/2)^2 - (j-N/2)^2));
            p = (i-M/2)/Depth(i,j);
            q = (j-N/2)/Depth(i,j);
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
            
        end
    end
end

%Source passed to algorithm used for reconstruction
s = [0,0];

figure;
imshow(mat2gray(Depth));

% Noisy image with gaussian noise
E_noise = imnoise(E,'gaussian',0,5);

save('DataFile1.mat', 'E', 's','radius','mask','boundary', 'p_init', 'q_init', 'f_init', 'g_init', 'E_noise', 'Depth');