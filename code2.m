% p,q estimate from single source

clear all

load('DataFile1.mat');

M = size(E,1);
N = size(E,2);

% Choose the En you want to use from E_noise(the one with gaussian noise)
% or E (noiseless) by commenting out the other one
% En = E_noise;
En = E;

% p_init, q_init contain the initial conditions (values on the boundary)
p_old = p_init;
q_old = q_init;

pn = zeros(size(En));
qn = zeros(size(En));

% Set the number of iterations, value of the regularization (lambda)
iters = 800;
lambda = 0.5  ;

% We use the iterative method derived from variational calculus to estimate the values of p,q
for kk = 1:iters,
    disp(kk)
    for i=2:(M-1),
        for j=2:(N-1),
            if(boundary(i,j)==0 && mask(i,j) ==1)
                % If the points lies not on the boundary, only then compute the value of p,q
                pn(i,j) = Ravg(p_old,i,j) + (1/lambda)*( En(i,j) - Rval(p_old(i,j), q_old(i,j), s))*Rp(p_old(i,j),q_old(i,j),s(1),s(2));
                qn(i,j) = Ravg(q_old,i,j) + (1/lambda)*( En(i,j) - Rval(p_old(i,j), q_old(i,j), s))*Rq(p_old(i,j),q_old(i,j),s(1),s(2));
            else 
                % Otherwise, since it is already on the boundary, we already have the correct value 
                pn(i,j) = p_old(i,j);
                qn(i,j) = q_old(i,j);
            end
        end
    end
    % Re assign the values for next iteration
    p_old = pn;
    q_old = qn;
end

save('DataFile2.mat', 'pn', 'qn');
