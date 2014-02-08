% Shape from Shading, EE 702, 2014
% Ashwin Kachhara, Sumeet Fefar

M=1000;
N=1000;

k=400;

Depth = zeros(M,N);
E = 0 * ones(M,N);
% p = 0 * ones(M,N);
% q = 0 * ones(M,N);
% ptemp = 0 * ones(M,N);
% qtemp = 0 * ones(M,N);
% f = 0 * ones(M,N);
% g = 0 * ones(M,N);
s=[0,0];

for i=1:M,
    for j=1:N,
        if(sqrt((i-M/2)^2 + (j-N/2)^2)<=k)
            Depth(i,j) = round(sqrt(k^2 - (i-M/2)^2 - (j-N/2)^2));
            p = (i-M/2)/Depth(i,j);
            q = (j-N/2)/Depth(i,j);
%             f(i,j) = (i-M/2)/(Depth(i,j)+ M/2);
%             g(i,j) = (j-N/2)/(Depth(i,j)+ N/2);
            temp = (s(1)*p + s(2)*q + 1)/sqrt( (s(2)^2 + s(1)^2 +1) * (p^2 + q^2 +1));
            if(temp>0)
                E(i,j) = temp;
            end
            
        end
    end
end

imshow(mat2gray(Depth));
imshow(mat2gray(E));

save('Esfile.mat', 'E', 's');



